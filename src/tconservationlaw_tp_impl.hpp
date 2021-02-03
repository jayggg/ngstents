#ifndef CONSERVATIONLAW_TP_IMPL
#define CONSERVATIONLAW_TP_IMPL

#include "conservationlaw.hpp"
#include "paralleldepend.hpp"
#include "tentsolver_impl.hpp"

template <typename EQUATION, int DIM, int COMP, int ECOMP>
void T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
CalcFluxTent (int tentnr, FlatMatrixFixWidth<COMP> u, FlatMatrixFixWidth<COMP> u0,
	      FlatMatrixFixWidth<COMP> flux, double tstar, LocalHeap & lh)
{
  static Timer tflux ("CalcFluxTent", 2);
  ThreadRegionTimer reg(tflux, TaskManager::GetThreadId());

  // note: tstar is not used
  const Tent & tent = tps->GetTent(tentnr);
  auto fedata = tent.fedata;
  if (!fedata) throw Exception("fedata not set");

  // these variables no longer exist
  //*(tent->time) = tent.timebot + tstar*(tent.ttop-tent.tbot);

  flux = 0.0;
  {
  for (int i : Range(tent.els))
    {
      HeapReset hr(lh);
      const DGFiniteElement<DIM> & fel =
	static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[i]);

      auto & simd_ir = *fedata->iri[i];
      auto & simd_mir = *fedata->miri[i];

      IntRange dn = tent.fedata->ranges[i];

      FlatMatrix<SIMD<double>> u_iptsa(COMP, simd_ir.Size(),lh);
      FlatMatrix<SIMD<double>> flux_iptsa(DIM*COMP, simd_ir.Size(),lh);
      FlatMatrix<SIMD<double>> flux_iptsa2(DIM*COMP, simd_ir.Size(),lh);

      fel.Evaluate (simd_ir, u.Rows(dn), u_iptsa);

      Cast().Flux(simd_mir, u_iptsa, flux_iptsa);

      FlatVector<SIMD<double>> di = fedata->adelta[i];
      for (auto k : Range(simd_ir.Size()))
        flux_iptsa.Col(k) *= simd_mir[k].GetWeight() * di(k);
      {
        for (size_t i = 0; i < COMP; i++)
          for (size_t j = 0; j < DIM; j++)
            flux_iptsa2.Row(i*DIM+j) = flux_iptsa.Row(j*COMP+i);

        fel.AddGradTrans (simd_mir, flux_iptsa2, flux.Rows(dn));
      }
    }
  }

  {
  for(int i : Range(tent.internal_facets))
    {
      HeapReset hr(lh);
      size_t elnr1 = fedata->felpos[i][0];
      size_t elnr2 = fedata->felpos[i][1];
      if(elnr2 != size_t(-1))
        {
          // inner facet
          const DGFiniteElement<DIM> & fel1 =
            static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[elnr1]);
          const DGFiniteElement<DIM> & fel2 =
            static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[elnr2]);

          IntRange dn1 = tent.fedata->ranges[elnr1];
          IntRange dn2 = tent.fedata->ranges[elnr2];

          auto & simd_ir_facet_vol1 = *fedata->firi[i][0];
          auto & simd_ir_facet_vol2 = *fedata->firi[i][1];

          int simd_nipt = simd_ir_facet_vol1.Size(); // IR's have the same size
          FlatMatrix<SIMD<double>> u1(COMP, simd_nipt, lh),
                                   u2(COMP, simd_nipt, lh);

          fel1.Evaluate(simd_ir_facet_vol1, u.Rows(dn1), u1);
          fel2.Evaluate(simd_ir_facet_vol2, u.Rows(dn2), u2);

          auto & simd_mir1 = *fedata->mfiri1[i];

          FlatMatrix<SIMD<double>> fn(COMP, simd_nipt, lh);
	  Cast().NumFlux (simd_mir1, u1, u2, fedata->anormals[i], fn);

          FlatVector<SIMD<double>> di = fedata->adelta_facet[i];
          for (size_t j : Range(simd_nipt))
            fn.Col(j) *= -1.0 * di(j) * simd_mir1[j].GetWeight();

          fel1.AddTrans(simd_ir_facet_vol1, fn, flux.Rows(dn1));
          fn *= -1.0;
          fel2.AddTrans(simd_ir_facet_vol2, fn, flux.Rows(dn2));
        }
      else
        {
          // boundary facet
          const DGFiniteElement<DIM> & fel1 =
            static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[elnr1]);

          IntRange dn1 = tent.fedata->ranges[elnr1];

          auto & simd_ir_facet_vol1 = *fedata->firi[i][0];

          int simd_nipt = simd_ir_facet_vol1.Size(); // IR's have the same size
          FlatMatrix<SIMD<double>> u1(COMP, simd_nipt, lh),
                                   u2(COMP, simd_nipt, lh);

          fel1.Evaluate(simd_ir_facet_vol1,u.Rows(dn1),u1);
          auto & simd_mir = *fedata->mfiri1[i];

          // check for boundary condition number
          int bc = bcnr[tent.internal_facets[i]];
          if (bc == 0) // outflow, use same values as on the inside
            {
              u2 = u1;
            }
          else if (bc == 1) // wall
            {
              Cast().u_reflect(simd_mir, u1, fedata->anormals[i], u2);
            }
          else if (bc == 2) // inflow, use initial data
            {
              fel1.Evaluate(simd_ir_facet_vol1, u0.Rows(dn1), u2);
            }
          else if (bc == 3) // transparent (wave)
            {
              Cast().u_transparent(simd_mir, u1, fedata->anormals[i], u2);
            }
          else
            throw Exception(string("no implementation for your") +
                  string(" chosen boundary condition number ") +
                  ToString(bc+1));

          FlatMatrix<SIMD<double>> fn(COMP, simd_nipt, lh);
	  Cast().NumFlux (simd_mir, u1, u2, fedata->anormals[i], fn);

          FlatVector<SIMD<double>> di = fedata->adelta_facet[i];
          for (size_t j : Range(simd_nipt))
            {
              auto fac = -1.0 * di(j) * simd_mir[j].GetWeight();
              fn.Col(j) *= fac;
            }

          fel1.AddTrans(simd_ir_facet_vol1, fn, flux.Rows(dn1));
        }
    }
  }

  for (int i : Range (tent.els))
    SolveM (tent, i, flux.Rows (tent.fedata->ranges[i]), lh);

}

template <typename EQUATION, int DIM, int COMP, int ECOMP>
void T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
CalcViscosityTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                   FlatMatrixFixWidth<COMP> ubnd, FlatVector<double> nu,
                   FlatMatrixFixWidth<COMP> visc, LocalHeap & lh)
{
  const Tent & tent = tps->GetTent(tentnr);

  // grad(u)*grad(v) - {du/dn} * [v] - {dv/dn} * [u] + alpha * p^2 / h * [u]*[v]
  double alpha = 4.0;
  alpha *= sqr(fes->GetOrder());

  auto fedata = tent.fedata;
  if (!fedata) throw Exception("fedata not set");

  visc = 0.0;
  for (size_t i : Range(tent.els))
    {
      HeapReset hr(lh);
      const DGFiniteElement<DIM> & fel =
	static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[i]);
      auto & simd_mir = *fedata->miri[i];

      IntRange dn = tent.fedata->ranges[i];

      FlatMatrix<SIMD<double>> gradu(DIM,simd_mir.Size(),lh);
      for(size_t j : Range(COMP))
        {
	  fel.EvaluateGrad(simd_mir,u.Col(j).Range(dn),gradu);
	  for(size_t k : Range(simd_mir.Size()))
	    gradu.Col(k) *= nu(i) * simd_mir[k].GetWeight();

	  fel.AddGradTrans(simd_mir,gradu,visc.Col(j).Range(dn));
        }

      //quick hack for other facets
      auto fnums = ma->GetElFacets (tent.els[i]);
      for(auto j : Range(fnums))
        if(!tent.internal_facets.Contains(fnums[j]))
          {
            int locfacetnr = j;
            auto & trafo = *fedata->trafoi[i];
            auto etfacet = ElementTopology::GetFacetType (trafo.GetElementType(),
                                                          locfacetnr);
            // note: need to add 1 to integration order
            SIMD_IntegrationRule fir (etfacet,2*fel.Order()+1);
            auto vnums = ma->GetElVertices (ElementId(VOL,tent.els[i]));
            Facet2ElementTrafo transform(trafo.GetElementType(), vnums);
            auto & simd_ir_facet_vol = transform(locfacetnr,fir,lh);
            auto & simd_mfir = trafo(simd_ir_facet_vol,lh);
            simd_mfir.ComputeNormalsAndMeasure(trafo.GetElementType(),locfacetnr);

            int simd_nipt = simd_ir_facet_vol.Size();
            FlatMatrix<SIMD<double>>
              u1(COMP, simd_nipt, lh), u2(COMP, simd_nipt,lh),
              jumpu(COMP, simd_nipt, lh), gradu(DIM, simd_nipt, lh),
              dudn(COMP, simd_nipt, lh), temp(DIM, simd_nipt, lh);

            fel.Evaluate(simd_ir_facet_vol,u.Rows(dn),u1);
            fel.Evaluate(simd_ir_facet_vol,ubnd.Rows(dn),u2);
            jumpu = u1-u2;

            FlatVector<SIMD<double>> fac(simd_nipt,lh);
            for(size_t k : Range(simd_nipt))
              {
                fac(k) = nu(i)*simd_mfir[k].GetWeight();
                jumpu.Col(k) *= fac(k);
              }

            auto normal = simd_mfir.GetNormals();
            for(size_t j : Range(COMP))
              {
                fel.EvaluateGrad(simd_mfir,u.Col(j).Range(dn),gradu);
                for(size_t k : Range(simd_nipt))
                  {
                    dudn(j,k) = -fac(k)*InnerProduct(gradu.Col(k),normal.Row(k));
                    temp.Col(k) = -jumpu(j,k) * normal.Row(k);
                  }
                fel.AddGradTrans(simd_mfir,temp,visc.Col(j).Range(dn));
              }
            auto h = fabs(simd_mfir[0].GetJacobiDet())/simd_mfir[0].GetMeasure();
            jumpu *= alpha / h;
            dudn += jumpu;
            fel.AddTrans(simd_ir_facet_vol,dudn,visc.Rows(dn));
          }
    }

  for(size_t i : Range(tent.internal_facets))
    {
      size_t elnr1 = fedata->felpos[i][0];
      size_t elnr2 = fedata->felpos[i][1];
      if(elnr2 != size_t(-1))
        {
          // inner facet
          HeapReset hr(lh);
          const DGFiniteElement<DIM> & fel1 =
            static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[elnr1]);
          const DGFiniteElement<DIM> & fel2 =
            static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[elnr2]);

          IntRange dn1 = tent.fedata->ranges[elnr1];
          IntRange dn2 = tent.fedata->ranges[elnr2];

          auto & simd_ir_facet_vol1 = *fedata->firi[i][0];
          auto & simd_ir_facet_vol2 = *fedata->firi[i][1];
          auto & simd_mir = *fedata->mfiri1[i];
          auto & simd_mir2 = *fedata->mfiri2[i];

          int simd_nipt = simd_ir_facet_vol1.Size(); // IR's have the same size
          FlatMatrix<SIMD<double>>
            u1(COMP, simd_nipt, lh), u2(COMP, simd_nipt, lh),
            gradu1(DIM, simd_nipt, lh), gradu2(DIM, simd_nipt, lh),
            jumpu(COMP, simd_nipt, lh), dudn(COMP, simd_nipt, lh),
            temp(DIM, simd_nipt, lh);

          fel1.Evaluate(simd_ir_facet_vol1,u.Rows(dn1),u1);
          fel2.Evaluate(simd_ir_facet_vol2,u.Rows(dn2),u2);
          jumpu = nu(elnr1)*u1-nu(elnr2)*u2;

          FlatVector<SIMD<double>> fac(simd_nipt,lh);
          for(size_t k : Range(simd_nipt))
            {
              fac(k) = simd_mir[k].GetWeight();
              jumpu.Col(k) *= fac(k);
            }

          auto normal = fedata->anormals[i];
          // *testout << normal << endl;
          for(size_t j : Range(COMP))
            {
              fel1.EvaluateGrad(simd_mir,u.Col(j).Range(dn1),gradu1);
	      fel2.EvaluateGrad(simd_mir2,u.Col(j).Range(dn2),gradu2);
	      for(size_t k : Range(simd_nipt))
		{
                  temp.Col(k) = -0.5 * jumpu(j,k) * normal.Col(k);
                  dudn(j,k) = InnerProduct(nu(elnr1)*gradu1.Col(k) +
                                           nu(elnr2)*gradu2.Col(k),normal.Col(k));
		  dudn(j,k) *= -0.5 * fac(k);
		}
	      fel1.AddGradTrans(simd_mir,temp,visc.Col(j).Range(dn1));
	      temp *= -1.0;
              fel2.AddGradTrans(simd_mir2,temp,visc.Col(j).Range(dn2));
            }
          auto h = fabs(simd_mir[0].GetJacobiDet())/simd_mir[0].GetMeasure();
          jumpu *= alpha / h;
          dudn += jumpu;
          fel1.AddTrans(simd_ir_facet_vol1,dudn,visc.Rows(dn1));
	  dudn *= -1.0;
	  fel2.AddTrans(simd_ir_facet_vol2,dudn,visc.Rows(dn2));
        }
    }

  for (int i : Range (tent.els))
    {
      SolveM (tent, i, fedata->adelta[i], visc.Rows (tent.fedata->ranges[i]), lh);
    }
}

template <typename EQUATION, int DIM, int COMP, int ECOMP>
void T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
CalcEntropyResidualTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                         FlatMatrixFixWidth<COMP> ut,
                         FlatMatrixFixWidth<ECOMP> res,
                         FlatMatrixFixWidth<COMP> u0, double tstar,
                         LocalHeap & lh)
{
  HeapReset hr(lh);

  const Tent & tent = tps->GetTent(tentnr);
  auto fedata = tent.fedata;
  if (!fedata) throw Exception("fedata not set");

  res = 0.0;
  for (int i : Range(tent.els))
    {
      HeapReset hr(lh);
      const DGFiniteElement<DIM> & fel =
	static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[i]);
      const SIMD_IntegrationRule & simd_ir = *fedata->iri[i];

      IntRange dn = tent.fedata->ranges[i];

      int ndof = fel.GetNDof();

      FlatVector<> uiel_coeff(ndof,lh);
      FlatVector<> cons_coeff(ndof,lh);

      FlatMatrix<SIMD<double>> ui(COMP,simd_ir.Size(),lh),
                               uti(COMP,simd_ir.Size(),lh);
      fel.Evaluate(simd_ir,u.Rows(dn),ui);
      fel.Evaluate(simd_ir,ut.Rows(dn),uti);

      FlatMatrix<AutoDiff<1,SIMD<double>>> adu(COMP,simd_ir.Size(),lh);
      for(size_t k : Range(COMP))
        for(size_t l : Range(simd_ir.Size()))
          {
            adu(k,l).Value() = ui(k,l);
            adu(k,l).DValue(0) = uti(k,l);
          }

      auto & simd_mir = *fedata->miri[i];

      auto gradbot = fedata->agradphi_bot[i];
      auto gradtop = fedata->agradphi_top[i];
      FlatMatrix<AutoDiff<1,SIMD<double>>> gradphi_mat(DIM, simd_mir.Size(), lh);
      for(size_t k : Range(DIM))
        for(size_t l : Range(simd_mir.Size()))
          {
            gradphi_mat(k,l).Value() = (1-tstar)*gradbot(k,l) + tstar*gradtop(k,l);
            gradphi_mat(k,l).DValue(0) = gradtop(k,l) - gradbot(k,l);
          }

      Cast().InverseMap(simd_mir, gradphi_mat, adu);
      FlatMatrix<SIMD<double>> Ei(ECOMP,simd_ir.Size(),lh),
                               Fi(DIM*ECOMP,simd_ir.Size(),lh);
      Cast().CalcEntropy(adu, gradphi_mat, Ei, Fi);

      FlatVector<SIMD<double>> di = fedata->adelta[i];
      for(size_t k : Range(simd_ir.Size()))
        {
          Ei(0,k) *= simd_mir[k].GetWeight();
          auto fac = -1.0 * simd_mir[k].GetWeight() * di(k);
          for(size_t l : Range(DIM))
            Fi(l,k) *= fac;
          if(ECOMP>1)
            throw Exception(
                "not yet implemented for more than one entropy function");
        }

      fel.AddTrans(simd_ir,Ei,res.Rows(dn));
      fel.AddGradTrans (simd_mir, Fi, res.Col(0).Range(dn));
    }

  FlatMatrixFixWidth<COMP> temp(u.Height(),lh);
  Cyl2Tent(tentnr, tstar, u, temp, lh);

  for(int i : Range(tent.internal_facets))
    {
      HeapReset hr(lh);
      size_t elnr1 = fedata->felpos[i][0];
      size_t elnr2 = fedata->felpos[i][1];

      if(elnr2 != size_t(-1))
        {
          // inner facet
          const DGFiniteElement<DIM> & fel1 =
            static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[elnr1]);
          const DGFiniteElement<DIM> & fel2 =
            static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[elnr2]);

          IntRange dn1 = tent.fedata->ranges[elnr1];
          IntRange dn2 = tent.fedata->ranges[elnr2];

          auto & simd_ir_facet_vol1 = *fedata->firi[i][0];
          auto & simd_ir_facet_vol2 = *fedata->firi[i][1];

          int simd_nipt = simd_ir_facet_vol1.Size(); // IR's have the same size
          FlatMatrix<SIMD<double>> u1(COMP, simd_nipt, lh),
                                   u2(COMP, simd_nipt, lh);

          fel1.Evaluate(simd_ir_facet_vol1,temp.Rows(dn1),u1);
          fel2.Evaluate(simd_ir_facet_vol2,temp.Rows(dn2),u2);

          FlatMatrix<SIMD<double>> Fn(ECOMP, simd_nipt, lh);
          Cast().EntropyFlux(u1,u2,fedata->anormals[i],Fn);

          auto & simd_mir = *fedata->mfiri1[i];
          FlatVector<SIMD<double>> di = fedata->adelta_facet[i];
          for (size_t j : Range(simd_nipt))
            {
              auto fac = di(j) * simd_mir[j].GetWeight();
              Fn.Col(j) *= fac;
            }

          fel1.AddTrans(simd_ir_facet_vol1,Fn,res.Rows(dn1));
          Fn *= -1.0;
          fel2.AddTrans(simd_ir_facet_vol2,Fn,res.Rows(dn2));
        }
      else
        {
          // boundary facet
          const DGFiniteElement<DIM> & fel1 =
            static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[elnr1]);

          IntRange dn1 = tent.fedata->ranges[elnr1];

          auto & simd_ir_facet_vol1 = *fedata->firi[i][0];

          int simd_nipt = simd_ir_facet_vol1.Size();
          FlatMatrix<SIMD<double>> u1(COMP, simd_nipt, lh),
                                   u2(COMP, simd_nipt, lh);

          fel1.Evaluate(simd_ir_facet_vol1,temp.Rows(dn1),u1);

	  auto & simd_mir1 = *fedata->mfiri1[i];
	  // set u2 dofs based on boundary condition number
          int bc = bcnr[tent.internal_facets[i]];
          if (bc == 0) // outflow, use same values as on the inside
            {
              u2 = u1;
            }
          else if (bc == 1) // wall
            {
              Cast().u_reflect(simd_mir1, u1, fedata->anormals[i], u2);
            }
          else if (bc == 2) // inflow, use initial data
            {
              fel1.Evaluate(simd_ir_facet_vol1, u0.Rows(dn1), u2);
            }
          else if (bc == 3) // transparent (wave)
            {
              Cast().u_transparent(simd_mir1, u1, fedata->anormals[i], u2);
            }
          else
            throw Exception(string(
               "no implementation for your chosen boundary condition number ")
                + ToString(bc+1));

          FlatMatrix<SIMD<double>> Fn(ECOMP, simd_nipt, lh);
          Cast().EntropyFlux(u1,u2,fedata->anormals[i],Fn);

	  FlatVector<SIMD<double>> di = fedata->adelta_facet[i];
          for (size_t j : Range(simd_nipt))
            {
              auto fac = di(j) * simd_mir1[j].GetWeight();
              Fn.Col(j) *= fac;
            }
          if(bc!=1)
            fel1.AddTrans(simd_ir_facet_vol1,Fn,res.Rows(dn1));
        }
    }
  for (int i : Range (tent.els))
    SolveM (tent, i, res.Rows (tent.fedata->ranges[i]), lh);
}

template <typename EQUATION, int DIM, int COMP, int ECOMP>
double T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
CalcViscosityCoefficientTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                              FlatMatrixFixWidth<ECOMP> res,
			      double tstar, LocalHeap & lh)
{
  const Tent & tent = tps->GetTent(tentnr);
  auto fedata = tent.fedata;
  if (!fedata) throw Exception("fedata not set");

  double nu_tent = 0.0;

  for (int i : Range(tent.els))
    {
      HeapReset hr(lh);
      ElementId ei (VOL, tent.els[i]);
      const DGFiniteElement<DIM> & fel =
	static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[i]);
      const SIMD_IntegrationRule & simd_ir = *fedata->iri[i];

      IntRange dn = tent.fedata->ranges[i];

      FlatMatrix<SIMD<double>> resi(ECOMP,simd_ir.Size(),lh);
      FlatMatrix<SIMD<double>> ui(COMP,simd_ir.Size(),lh);

      auto & simd_mir = *fedata->miri[i];
      double hi = pow(simd_mir[0].GetMeasure()[0]/DIM,1.0/DIM);
      if( fel.Order() > 0 ) hi /= fel.Order();

      fel.Evaluate(simd_ir,u.Rows(dn),ui);
      fel.Evaluate(simd_ir,res.Rows(dn),resi);

      FlatVector<SIMD<double>> di = fedata->adelta[i];
      for(size_t k : Range(simd_ir.Size()))
        resi.Col(k) /= di(k);

      // clear overhead
      FlatMatrix<double> temp(ECOMP,simd_ir.Size()*SIMD<double>::Size(),
                              &resi(0,0)[0]);
      temp.Cols(simd_ir.GetNIP(),simd_ir.Size()*SIMD<double>::Size()) = 0.0;
      FlatMatrix<double> tempu(COMP,simd_ir.Size()*SIMD<double>::Size(),
                               &ui(0,0)[0]);
      tempu.Cols(simd_ir.GetNIP(),simd_ir.Size()*SIMD<double>::Size()) = 0.0;

      FlatMatrix<SIMD<double>> gradphi_mat(DIM, simd_mir.Size(), lh);
      gradphi_mat = (1-tstar)*fedata->agradphi_bot[i] +
                    tstar*fedata->agradphi_top[i];

      Cast().InverseMap(simd_mir, gradphi_mat, ui);
      Cast().CalcViscCoeffEl(simd_mir, ui, resi, hi, nu(ei.Nr()));

      if(nu(ei.Nr()) > nu_tent)
	nu_tent = nu(ei.Nr());
    }
  return nu_tent;
}

////////////////////////////////////////////////////////////////
// implementations of maps 
////////////////////////////////////////////////////////////////

template <typename EQUATION, int DIM, int COMP, int ECOMP>
void T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
Cyl2Tent (int tentnr, double tstar,
	  FlatMatrixFixWidth<COMP> uhat, FlatMatrixFixWidth<COMP> u,
	  LocalHeap & lh)
{
  static Timer tcyl2tent ("Cyl2Tent", 2);
  ThreadRegionTimer reg(tcyl2tent, TaskManager::GetThreadId());

  const Tent & tent = tps->GetTent(tentnr);

  auto fedata = tent.fedata;
  if (!fedata) throw Exception("fedata not set");

  for (size_t i : Range(tent.els))
    {
      HeapReset hr(lh);
      const DGFiniteElement<DIM> & fel =
        static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[i]);

      auto & simd_mir = *fedata->miri[i];
      IntRange dn = tent.fedata->ranges[i];

      FlatMatrix<SIMD<double>> u_ipts(COMP, simd_mir.Size(),lh);
      for(size_t k = 0; k < COMP; k++)
        fel.Evaluate (simd_mir.IR(), uhat.Col(k).Range(dn), u_ipts.Row(k));

      FlatMatrix<SIMD<double>> gradphi_mat(DIM, simd_mir.Size(), lh);
      gradphi_mat = (1-tstar)*fedata->agradphi_bot[i] +
                    tstar*fedata->agradphi_top[i];

      Cast().InverseMap(simd_mir, gradphi_mat, u_ipts);

      for(size_t k = 0; k < simd_mir.Size(); k++)
        u_ipts.Col(k) *= simd_mir[k].GetWeight();

      u.Rows(dn) = 0.0;
      for(size_t k = 0; k < COMP; k++)
        fel.AddTrans (simd_mir.IR(), u_ipts.Row(k), u.Col(k).Range(dn));

      SolveM (tent, i, u.Rows (dn), lh);
    }
}

template <typename EQUATION, int DIM, int COMP, int ECOMP>
void T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
ApplyM1 (int tentnr, double tstar, FlatMatrixFixWidth<COMP> u,
         FlatMatrixFixWidth<COMP> res, LocalHeap & lh)
{
  static Timer tapplym1 ("ApplyM1", 2);
  ThreadRegionTimer reg(tapplym1, TaskManager::GetThreadId());

  const Tent & tent = tps->GetTent(tentnr);

  auto fedata = tent.fedata;
  if (!fedata) throw Exception("fedata not set");

  res = 0.0;
  for (int i : Range(tent.els))
    {
      HeapReset hr(lh);
      const DGFiniteElement<DIM> & fel =
	static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[i]);

      const SIMD_IntegrationRule & simd_ir = *fedata->iri[i];
      IntRange dn = tent.fedata->ranges[i];

      FlatMatrix<SIMD<double>> u_ipts(COMP, simd_ir.Size(), lh),
                               temp(COMP, simd_ir.Size(), lh);
      FlatMatrix<SIMD<double>> flux(COMP*DIM, simd_ir.Size(), lh);
      FlatMatrix<SIMD<double>> graddelta_mat(DIM, simd_ir.Size(), lh);
      graddelta_mat = fedata->agradphi_top[i] - fedata->agradphi_bot[i];

      fel.Evaluate (simd_ir, u.Rows(dn), u_ipts);

      auto & simd_mir = *fedata->miri[i];

      Cast().Flux(simd_mir,u_ipts,flux);
      for(size_t j : Range(simd_ir.Size()))
        for(size_t l : Range(COMP))
          {
            SIMD<double> hsum(0.0);
            for(size_t k : Range(DIM))
              {
                auto graddelta = graddelta_mat(k,j) * simd_mir[j].GetWeight();
                hsum += graddelta * flux(COMP*k+l,j);
              }
            temp(l,j) = hsum;
          }
      fel.AddTrans(simd_ir,temp,res.Rows(dn));

      SolveM (tent, i, res.Rows (dn), lh);
    }
}

template <typename EQUATION, int DIM, int COMP, int ECOMP>
void T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
Tent2Cyl (int tentnr, double tstar,
	  FlatMatrixFixWidth<COMP> u, FlatMatrixFixWidth<COMP> uhat,
          bool solvemass, LocalHeap & lh)
{
  static Timer ttent2cyl ("Tent2Cyl", 2);
  ThreadRegionTimer reg(ttent2cyl, TaskManager::GetThreadId());

  const Tent & tent = tps->GetTent(tentnr);

  auto fedata = tent.fedata;
  if (!fedata) throw Exception("fedata not set");

  uhat = 0.0;
  for (int i : Range(tent.els))
    {
      HeapReset hr(lh);
      const DGFiniteElement<DIM> & fel =
	static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[i]);
      const SIMD_IntegrationRule & simd_ir = *fedata->iri[i];

      IntRange dn = tent.fedata->ranges[i];

      FlatMatrix<SIMD<double>> u_ipts(COMP, simd_ir.Size(), lh);
      FlatMatrix<SIMD<double>> flux(COMP*DIM, simd_ir.Size(), lh),
                               res(COMP, simd_ir.Size(), lh);
      FlatMatrix<SIMD<double>> gradphi_mat(DIM, simd_ir.Size(), lh);
      gradphi_mat = (1-tstar)*fedata->agradphi_bot[i] +
                    tstar*fedata->agradphi_top[i];

      fel.Evaluate (simd_ir, u.Rows(dn), u_ipts);

      auto & simd_mir = *fedata->miri[i];
      Cast().Flux(simd_mir,u_ipts,flux);

      for(size_t j : Range(simd_ir.Size()))
        for(size_t l : Range(COMP))
          {
            SIMD<double> hsum(0.0);
            for(size_t k : Range(DIM))
              {
                auto gradphi = gradphi_mat(k,j) * simd_mir[j].GetWeight();
                hsum += gradphi * flux(COMP*k+l,j);
              }
            res(l,j) = u_ipts(l,j) * simd_mir[j].GetWeight() - hsum;
          }

      fel.AddTrans(simd_ir,res,uhat.Rows(dn));
      if(solvemass)
	SolveM (tent, i, uhat.Rows (dn), lh);
    }
}

////////////////////////////////////////////////////////////////
// time stepping methods 
////////////////////////////////////////////////////////////////

template <typename EQUATION, int DIM, int COMP, int ECOMP>
void T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
Propagate(LocalHeap & lh)
{
  static Timer tprop ("Propagate new", 2); RegionTimer reg(tprop);
  RunParallelDependency
    (tent_dependency, [&] (int i)
     {
       LocalHeap slh = lh.Split();  // split to threads
       // const Tent & tent = tps->GetTent(i);
       tentsolver->PropagateTent(i, *u, *uinit, slh);
     });
}

template <typename EQUATION, int DIM, int COMP, int ECOMP>
void T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
PropagateSAT(int stages, int substeps,
	     BaseVector & hu, BaseVector & hu0,
	     LocalHeap & lh)
{
  static Timer t_SAT ("PropagateSAT", 2); RegionTimer reg(t_SAT);
  static Timer t_SATtent ("PropagateSAT - tent", 2);

  RunParallelDependency 
    (tent_dependency, [&] (int i)
     {
       size_t tid = TaskManager::GetThreadId();
       ThreadRegionTimer reg(t_SATtent, tid);
       RegionTracer reg1(tid, t_SATtent, i);

       LocalHeap slh = lh.Split();  // split to threads

       const Tent & tent = tps->GetTent(i);
       tent.fedata = new (slh) TentDataFE(tent, *fes, *ma, slh);

       int ndof = tent.fedata->nd;
       FlatMatrixFixWidth<COMP> local_uhat(ndof,slh);
       FlatMatrixFixWidth<COMP> local_u0(ndof,slh);
       FlatMatrixFixWidth<COMP> local_u0temp(ndof,slh);
       hu.GetIndirect(tent.fedata->dofs, AsFV(local_uhat));
       hu0.GetIndirect(tent.fedata->dofs, AsFV(local_u0));
       local_u0temp = local_u0;

       FlatMatrixFixWidth<COMP> local_uhat1(ndof,slh);
       FlatMatrixFixWidth<COMP> local_u(ndof,slh);
       FlatMatrixFixWidth<COMP> local_help(ndof,slh);

       double taustar = 1.0/substeps;
       for (int j = 0; j < substeps; j++)
	 {
	   local_uhat1 = local_uhat;
	   double fac = 1.0;
	   local_u0 = local_u0temp;
	   for(int k : Range(1,stages+1))
	     {
	       Cyl2Tent(i, j*taustar, local_uhat1, local_u, slh);
	       CalcFluxTent(i, local_u, local_u0, local_uhat1, j*taustar, slh);
               local_uhat1 *= 1.0/k;
	       fac *= taustar;
	       local_uhat += fac*local_uhat1;           

	       if(k < stages)
	   	 {
		   ApplyM1(i, j*taustar, local_u, local_help, slh);
	   	   local_uhat1 += local_help;
	   	 }
	       local_u0 = 0.0;
	     }
         }
       hu.SetIndirect(tent.fedata->dofs, AsFV(local_uhat));
       tent.fedata = nullptr;
     });
}

template <typename EQUATION, int DIM, int COMP, int ECOMP>
void T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>::
PropagateSARK(int stages, int substeps, BaseVector & hu, BaseVector & hu_init, LocalHeap & lh)
{
  shared_ptr<BaseVector> hres = (ECOMP > 0) ? gfres->GetVectorPtr() : nullptr;
  
  RunParallelDependency
    (tent_dependency, [&] (int i)
     {
       *testout << "tent " << i << endl;
       LocalHeap slh = lh.Split();  // split to threads
       const Tent & tent = tps->GetTent(i);
       tent.fedata = new (slh) TentDataFE(tent, *fes, *ma, slh);

       int ndof = (tent.fedata)->nd;
       double taustar = 1.0/substeps;

       // tent.fedata = new (slh) TentDataFE(tent, *fes, *ma, slh);
       // tent->InitTent(gftau);

       FlatMatrixFixWidth<COMP> local_u0(ndof,slh);
       FlatMatrixFixWidth<COMP> local_Gu0(ndof,slh);
       FlatMatrixFixWidth<COMP> local_init(ndof,slh);

       hu.GetIndirect(tent.fedata->dofs, AsFV(local_Gu0));
       hu_init.GetIndirect(tent.fedata->dofs, AsFV(local_init));

       FlatMatrixFixWidth<COMP> local_u1(ndof,slh);
       FlatMatrixFixWidth<COMP> local_u1_half(ndof,slh);
       FlatMatrixFixWidth<COMP> local_u(ndof,slh);
       FlatMatrixFixWidth<COMP> local_help(ndof,slh);
       FlatMatrixFixWidth<COMP> local_flux(ndof,slh), local_flux1(ndof,slh);

       FlatMatrixFixWidth<COMP> U0(ndof,slh);
       FlatMatrixFixWidth<COMP> U1(ndof,slh);
       FlatMatrixFixWidth<COMP> U2(ndof,slh);
       FlatMatrixFixWidth<COMP> U3(ndof,slh);
       FlatMatrixFixWidth<COMP> u0(ndof,slh);
       FlatMatrixFixWidth<COMP> u1(ndof,slh);
       FlatMatrixFixWidth<COMP> u2(ndof,slh);
       FlatMatrixFixWidth<COMP> u3(ndof,slh);
       FlatMatrixFixWidth<COMP> M1u0(ndof,slh);
       FlatMatrixFixWidth<COMP> M1u1(ndof,slh);
       FlatMatrixFixWidth<COMP> M1u2(ndof,slh);
       FlatMatrixFixWidth<COMP> M1u3(ndof,slh);
       FlatMatrixFixWidth<COMP> fu0(ndof,slh);
       FlatMatrixFixWidth<COMP> fu1(ndof,slh);
       FlatMatrixFixWidth<COMP> fu2(ndof,slh);
       FlatMatrixFixWidth<COMP> fu3(ndof,slh);
       FlatMatrixFixWidth<COMP> Uhat(ndof,slh);

       FlatMatrixFixWidth<COMP> dUhatdt(ndof,slh);
       FlatMatrixFixWidth<ECOMP> res(ndof,slh);
       FlatVector<> local_nu(tent.els.Size(),slh);

       double tau_tent = tent.ttop - tent.tbot;
       double h_tent = 0;
       for (int j : Range(tent.els))
         h_tent = max2(h_tent, tent.fedata->mesh_size[j]);

       int order = max(1,fes->GetOrder());

       double tau_visc1 = sqr (h_tent / sqr(order));

       // calc |u|_M0 on advancing front
       Cyl2Tent (i, 0, local_Gu0, local_u0, slh);
       Tent2Cyl(i, 0, local_u0, local_help, false, slh);

       double norm_bot = InnerProduct(AsFV(local_u0),AsFV(local_help));

       //TODO: implement SARK for different number to stages
       for (int j = 0; j < substeps; j++)
	 {
           // third order
           U0 = local_Gu0;
           Cyl2Tent (i, j*taustar, U0, u0, slh);
           ApplyM1(i, j*taustar, u0, M1u0, slh);
           CalcFluxTent(i, u0, local_init, fu0, j*taustar, slh);
           U1 = U0 + 0.5*taustar*(M1u0+fu0);

           Cyl2Tent (i, j*taustar, U1, u1, slh);
           ApplyM1(i, j*taustar, u1, M1u1, slh);
           CalcFluxTent(i, u1, local_init, fu1, (j+0.5)*taustar, slh);
           U2 = U0 + taustar*(4*M1u1-3*M1u0+2*fu1-fu0);

           Cyl2Tent (i, j*taustar, U2, u2, slh);
           CalcFluxTent(i, u2, local_init, fu2, (j+1)*taustar, slh);

           Uhat = U0 + taustar * (1.0/6.0 * fu0 + 2.0/3.0 * fu1 + 1.0/6.0 * fu2);
           local_Gu0 = Uhat;
           // for viscosity
           dUhatdt = fu0;
           if (ECOMP > 0)
	     {
               /////// use dUhatdt as approximation at the final time
	       // U0 = local_Gu0;
               // Cyl2Tent (i, (j+1)*taustar, local_Gu0, local_help, slh);
               // CalcFluxTent(i, local_help, local_init, dUhatdt,
               //              (j+1)*taustar, slh);
	       // CalcEntropyResidualTent(i, U0, dUhatdt, res, local_init,
               //                         (j+1)*taustar, slh);
               // hres->SetIndirect(tent.dofs,AsFV(res));
	       // double nu_tent = CalcViscosityCoefficientTent(
               //                        i, U0, res, (j+1)*taustar, slh);
	       /////// use dUhatdt as approximation at the initial time
               CalcEntropyResidualTent(i, U0, dUhatdt, res, local_init,
                                       j*taustar, slh);
               hres->SetIndirect(tent.fedata->dofs,AsFV(res));
	       double nu_tent = CalcViscosityCoefficientTent(i, U0, res,
                                                             j*taustar, slh);

               local_nu = nu_tent;
	       double steps_visc = (40*tau_tent*nu_tent/tau_visc1)/substeps;
               if (steps_visc > 0.2)
                 {
                   steps_visc = max(1.0,ceil(steps_visc));
                   double tau_visc = taustar/steps_visc;
		   // store boundary conditions in local_help
		   Cyl2Tent (i, (j+1)*taustar, local_Gu0, local_u, slh);
		   local_help = local_u;
                   for (int k = 0; k < steps_visc; k++)
                     {
		       CalcViscosityTent (i, local_u, local_help, local_nu,
                                          local_flux, slh);
                       local_u -= tau_visc * local_flux;
                     }
                   Tent2Cyl(i, (j+1)*taustar, local_u, local_Gu0, true, slh);
                 }
	     }
         }

       // calc |u|_M1 norm on advancing front
       Cyl2Tent (i, 1, local_Gu0, local_u0, slh);
       Tent2Cyl(i, 1, local_u0, local_help, false, slh);
       double norm_top = InnerProduct(AsFV(local_u0),AsFV(local_help));
       if(norm_bot>0)
         *testout << "top/bot = " << norm_top/norm_bot << endl;
       else
         *testout << "bot, top = " << norm_bot << ", " << norm_top << endl;
       if(norm_top/norm_bot > 1.0)
         *testout << "bot, top : " << norm_bot << ", " << norm_top << endl;
       if(norm_top/norm_bot < 0.9)
         *testout << "bot, top : " << norm_bot << ", " << norm_top << endl;

       hu.SetIndirect(tent.fedata->dofs, AsFV(local_Gu0));
       tent.fedata = nullptr;
       // Tent no longer has this method
       // tps.GetTent(i)->SetFinalTime();
     });
}
#endif // CONSERVATIONLAW_TP_IMPL
