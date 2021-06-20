#include "vis3d.hpp"

void Visualization3D::SetInitialHd(
    shared_ptr<GridFunction> gfu,
    shared_ptr<GridFunction> hdgf, LocalHeap & lh)
{
    // Assuming no Complex valued spaces for now
    auto fes = gfu->GetFESpace();
    auto hdfes = hdgf->GetFESpace();
    int order = hdfes->GetOrder();

    shared_ptr<MeshAccess> ma = fes->GetMeshAccess();

    // Note that tmph1 and vtmp are attributes of the class
    if  (ma->GetNPeriodicIdentifications() > 0) {

      auto tmph1a = CreateFESpace("h1ho", ma, Flags()
                            .SetFlag("order", order)
                            .SetFlag("dim", fes->GetDimension()));
      tmph1a->Update();
      tmph1a->FinalizeUpdate();
      Flags flags = tmph1a->GetFlags();
      Array<int> emptyarray;
      shared_ptr<Array<int>> used_idnrs = make_shared<Array<int>>(emptyarray);
      tmph1 = make_shared<PeriodicFESpace>(tmph1a, flags, used_idnrs);
    }
    else
      tmph1 = CreateFESpace("h1ho", ma, Flags()
                            .SetFlag("order", order)
                            .SetFlag("dim", fes->GetDimension()));

    tmph1->Update();
    tmph1->FinalizeUpdate();
    auto gftmp = CreateGridFunction(tmph1,"gftmp",Flags().SetFlag("novisual"));
    gftmp->Update();
    vtmp = gftmp->GetVectorPtr(0);

    ngcomp::SetValues (gfu, *gftmp, VOL, NULL, lh, false, true, 0);

    // Transfer dof values from vtmp to vhd.
    shared_ptr<BaseVector> vhd = hdgf->GetVectorPtr(0);
    int entrysize = vtmp->EntrySize();
    for (auto i : IntRange(vtmp->Size())) {
      int idx = (*(idx3d[0]))[i][0];
      if(entrysize == 1)
      {
         auto fv = vtmp->FVDouble();
         vhd->Range(idx,idx+1) = fv[i];
      }
      else
      {
         auto sv = vtmp->SV<double>();
         vhd->SV<double>()(idx) = sv(i);
      }
    }
}


void Visualization3D::SetForTent(
    Tent &tent, shared_ptr<GridFunction> gfu,
    shared_ptr<GridFunction> hdgf, LocalHeap & lh)
{
    auto fes = gfu->GetFESpace();
    auto hdfes = hdgf->GetFESpace();
    int order = hdfes->GetOrder();
    shared_ptr<MeshAccess> ma = fes->GetMeshAccess();
    int NV = ma->GetNV();
    int dim = tmph1->GetDimension();

    Array<int> cnti(tmph1->GetNDof());
    cnti = 0;

    DifferentialOperator * diffop = tmph1->GetEvaluator(VOL).get();
    shared_ptr<BilinearFormIntegrator> bli = tmph1->GetIntegrator(VOL);
    shared_ptr<BilinearFormIntegrator> single_bli = bli;
    if (dynamic_pointer_cast<BlockBilinearFormIntegrator> (single_bli))
      single_bli = dynamic_pointer_cast<BlockBilinearFormIntegrator> (single_bli)->BlockPtr();

    // So we can restore this afterwards
    bool bli_uses_simd = bli->SimdEvaluate(), sbli_uses_simd = single_bli->SimdEvaluate();

    int dimflux = diffop ? diffop->Dim() : bli->DimFlux();

    if (gfu -> Dimension() != dimflux)
      throw Exception(string("Error in SetValues: gridfunction-dim = ") + ToString(dimflux) +
                      ", but coefficient-dim = " + ToString(gfu->Dimension()));

    auto cachecfs = FindCacheCF (*gfu);

    // zeroing the vector vtmp would create a race condition. Instead
    // zero out only vector components associated with our tent's elements.
    IterateElements
      (*tmph1, VOL, lh,
       [&] (FESpace::Element ei, LocalHeap & lh)
       {
         if (tent.els.Contains(ei.Nr()))
         {
            const FiniteElement & fel = tmph1->GetFE (ei, lh);
            FlatVector<> elzero(fel.GetNDof() * dim, lh);
            elzero = 0.0;
            vtmp->SetIndirect (ei.GetDofs(), elzero);
         }
       });
    // TODO: Is there some way to iterate over only tent elements?
    IterateElements
      (*tmph1, VOL, lh,
       [&] (FESpace::Element ei, LocalHeap & lh)
       {
         if (tent.els.Contains(ei.Nr()))
         {
            const FiniteElement & fel = tmph1->GetFE (ei, lh);
            int ndof = fel.GetNDof();
            const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

            FlatVector<> elflux(ndof * dim, lh);
            FlatVector<> elfluxi(ndof * dim, lh);
            FlatVector<> fluxi(dimflux, lh);

            SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
            FlatMatrix<SIMD<double>> mfluxi(dimflux, ir.Size(), lh);

            auto & mir = eltrans(ir, lh);

            ProxyUserData ud;
            const_cast<ElementTransformation&>(eltrans).userdata = &ud;
            PrecomputeCacheCF (cachecfs, mir, lh);

            gfu->Evaluate (mir, mfluxi);

            for (size_t j : Range(ir))
              mfluxi.Col(j) *= mir[j].GetWeight();

            elflux = 0.0;
            diffop -> AddTrans (fel, mir, mfluxi, elflux);

            if (dim > 1)
              {
                FlatMatrix<double> elmat(fel.GetNDof(), lh);
                single_bli->CalcElementMatrix (fel, eltrans, elmat, lh);
                FlatCholeskyFactors<double> invelmat(elmat, lh);

                for (int j = 0; j < dim; j++)
                  invelmat.Mult (elflux.Slice (j,dim), elfluxi.Slice (j,dim));
              }
            else
              {
                FlatMatrix<double> elmat(fel.GetNDof(), lh);
                bli->CalcElementMatrix (fel, eltrans, elmat, lh);

                // Transform solution inverse instead
                CalcLDL<double,ColMajor> (Trans(elmat));
                elfluxi = elflux;
                SolveLDL<double,ColMajor> (Trans(elmat), elfluxi);
              }

            tmph1->TransformVec (ei, elfluxi, TRANSFORM_SOL_INVERSE);

            vtmp->GetIndirect (ei.GetDofs(), elflux);
            elfluxi += elflux;
            vtmp->SetIndirect (ei.GetDofs(), elfluxi);

            for (auto d : ei.GetDofs())
              if (IsRegularDof(d)) cnti[d]++;
         }
    });

    // Restore simd evaluate
    bli->SetSimdEvaluate(bli_uses_simd);
    single_bli->SetSimdEvaluate(sbli_uses_simd);

    VectorMem<10,double> fluxi(dim);
    ArrayMem<int,1> dnums(1);

    ParallelForRange
      (cnti.Size(), [&] (IntRange r)
       {
         VectorMem<10,double> fluxi(dim);
         ArrayMem<int,1> dnums(1);
         for (auto i : r)
           if (cnti[i])
             {
               dnums[0] = i;
               vtmp->GetIndirect (dnums, fluxi);
               fluxi /= cnti[i];
               vtmp->SetIndirect (dnums, fluxi);
             }
       });

    // Transfer dof values from vtmp to vhd
    Array<int> vtmp_nrs;
    Array<int> vhd_nrs;
    vtmp_nrs.Append(tent.vertex);

    if (order > 1) {
      // edge dofs offset by NV
      for (auto f : tent.internal_facets)
        vtmp_nrs.Append(NV+f);
    }

    for (auto i : vtmp_nrs)
      vhd_nrs.Append((*(idx3d[tent.level+1]))[i][0]);

    shared_ptr<BaseVector> vhd = hdgf->GetVectorPtr(0);
    int entrysize = vtmp->EntrySize();

    if(entrysize == 1)
    {
      for (auto i : IntRange(vtmp_nrs.Size()))
      {
        auto fv = vtmp->FVDouble();
        auto vi = vhd_nrs[i];
        vhd->Range(vi,vi+1) = fv[vtmp_nrs[i]];
      }
    }
    else
    {
      for (auto i : IntRange(vtmp_nrs.Size()))
      {
        auto sv = vtmp->SV<double>();
        vhd->SV<double>()(vhd_nrs[i]) = sv(vtmp_nrs[i]);
      }
    }
}
