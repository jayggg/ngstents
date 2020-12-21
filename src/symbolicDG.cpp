#include <comp.hpp>
#include <python_ngstd.hpp>


using namespace ngcomp;


static LocalHeap lh(100000, "applyDG", true);


class MTP_UserData : public ProxyUserData
{
public:
  FlatVector<> gradphi;
  using ProxyUserData::ProxyUserData;
};


void ApplyDGOperator (FESpace & fes,
                      CoefficientFunction & flux,
                      CoefficientFunction & numflux,
                      BaseVector & x, BaseVector & y)
{
  static Timer t("ApplyDGOp");
  static Timer tvol("ApplyDGOpVOL");
  static Timer tfac("ApplyDGOpFAC");
  static Timer tfac_apply("ApplyDGOpFAC - apply");
  static Timer tfac_applyt("ApplyDGOpFAC - applytran");
  static Timer tfac_getfe("ApplyDGOpFAC - getFE");

  RegionTimer reg(t);

  auto ma = fes.GetMeshAccess();

  Array<ProxyFunction*> trial_proxies, numflux_trial_proxies;

  flux.TraverseTree
    ( [&] (CoefficientFunction & nodecf)
      {
	// cout << "node-type = " << typeid(nodecf).name() << endl;
        auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
        if (proxy)
          {
            if (!proxy->IsTestFunction())
              if (!trial_proxies.Contains(proxy))
                trial_proxies.Append (proxy);
          }
      });

  numflux.TraverseTree
    ( [&] (CoefficientFunction & nodecf)
      {
        // cout << "node-type = " << typeid(nodecf).name() << endl;
        auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
        if (proxy)
          {
            if (!proxy->IsTestFunction())
              if (!numflux_trial_proxies.Contains(proxy))
                numflux_trial_proxies.Append (proxy);
          }
      });

  int num_equ = flux.Dimensions()[0];
  int dim_space = ma->GetDimension();

  // cout << "trial_proxies: " << trial_proxies << endl;
  // cout << "numflux_trial_proxies: " << numflux_trial_proxies << endl;
  // cout << "num_equ = " << num_equ << ", dim_space = " << dim_space << endl;
  // cout << "flux.Dimensions = " << flux.Dimensions() << endl;

  if (flux.Dimensions().Size() != 2)
    throw Exception ("flux should be a matrix");
  if (flux.Dimensions()[1] != dim_space)
    throw Exception ("what are you doing ???");
  if (fes.GetDimension() != num_equ)
    throw Exception ("space dimensions and flux rows don't match");

  shared_ptr<DifferentialOperator> gradient = fes.GetFluxEvaluator(VOL);
  // needed for trace
  shared_ptr<DifferentialOperator> evaluator = fes.GetEvaluator(VOL);

  y = 0.0;

  tvol.Start();

  IterateElements    // element applications in parallel
    (fes, VOL, lh,
     [&] (FESpace::Element el, LocalHeap & lh)
     {
       const FiniteElement & fel = el.GetFE();
       FlatArray<int> dnums = el.GetDofs();

       FlatVector<> elx(dnums.Size()*fes.GetDimension(), lh);
       FlatVector<> ely(dnums.Size()*fes.GetDimension(), lh);

       x.GetIndirect (dnums, elx);

       SIMD_IntegrationRule ir(fel.ElementType(), 2*fel.Order());
       const ElementTransformation & trafo = el.GetTrafo();
       auto & mir = trafo(ir, lh);

       MTP_UserData ud(trial_proxies.Size(), lh);
       const_cast<ElementTransformation&>(trafo).userdata = &ud;
       ud.fel = &fel;
       // ud.elx = &elx;
       // ud.lh = &lh;

       Vec<2> dummy_gradphi(1,0);
       ud.gradphi.AssignMemory (2, &dummy_gradphi(0));

       for (ProxyFunction * proxy : trial_proxies)
         ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
       for (ProxyFunction * proxy : trial_proxies)
         proxy->Evaluator()->Apply(fel, mir, elx, ud.GetAMemory(proxy));

       ely = 0;

       FlatMatrix<SIMD<double>> flux_pnt_values(flux.Dimension(), ir.Size(), lh);
       flux.Evaluate (mir, flux_pnt_values);

       for (size_t i = 0; i < flux_pnt_values.Height(); i++)
         {
           auto row = flux_pnt_values.Row(i);
           for (size_t j = 0; j < row.Size(); j++)
             row(j) *= mir[j].GetWeight();
         }

       gradient->AddTrans(fel, mir, flux_pnt_values, ely);
       y.AddIndirect (dnums, ely);
     });

  tvol.Stop();

  tfac.Start();

  for (auto colfacets : fes.FacetColoring())
    ParallelForRange
      (colfacets.Size(), [&] (IntRange r)
       {
         LocalHeap slh = lh.Split();
         LocalHeap & lh = slh;
         Array<int> elnums(2, lh), elnums_per(2, lh), // fnumsm(6, lh),
	   fnumsp(6, lh), vnumsm(8, lh), vnumsp(8, lh);

         for (size_t i : r)
           {
             HeapReset hr(slh);

             int facet = colfacets[i];
             ma->GetFacetElements (facet, elnums);
             if (elnums.Size() == 0) continue; // coarse facets

             int elm = elnums[0];
	     // ma->GetElFacets (elm,fnumsm);
	     auto fnumsm = ma->GetElFacets(elm);

             int facnrm = fnumsm.Pos(facet);

             ElementId eim(VOL, elm);

             if (elnums.Size() < 2)
               {
                 // boundary terms tbd
                 /*
                   ma->GetFacetSurfaceElements (facet, elnums);
                   int sel = elnums[0];
                   ElementId sei(BND, sel);

                   const FiniteElement & fel = fespace->GetFE (ei1, lh);
                   Array<int> dnums(fel.GetNDof(), lh);
                   ma->GetElVertices (el1, vnums1);
                   ma->GetSElVertices (sel, vnums2);

                   ElementTransformation & eltrans = ma->GetTrafo (ei1, lh);
                   ElementTransformation & seltrans = ma->GetTrafo (sei, lh);

                   fespace->GetDofNrs (ei1, dnums);

                   for (int j = 0; j < NumIntegrators(); j++)
                   {
                   const BilinearFormIntegrator & bfi = *parts[j];

                   if (!bfi.BoundaryForm()) continue;
                   if (!bfi.SkeletonForm()) continue;
                   if (bfi.GetDGFormulation().element_boundary) continue;
                   if (!bfi.DefinedOn (seltrans.GetElementIndex())) continue;

                   FlatVector<SCAL>
                     elx(dnums.Size()*this->fespace->GetDimension(), lh),
                     ely(dnums.Size()*this->fespace->GetDimension(), lh);
                   x.GetIndirect(dnums, elx);

                   dynamic_cast<const FacetBilinearFormIntegrator&>(bfi).
                   ApplyFacetMatrix (fel,facnr1,eltrans,vnums1, seltrans,
                                     vnums2, elx, ely, lh);
                   y.AddIndirect(dnums, ely);
                   } //end for (numintegrators)

                   continue;
                 */
                 continue;
               } // end if boundary facet


             int elp = elnums[1];
             ElementId eip(VOL, elp);

             auto fnumsp = ma->GetElFacets(elp);
             int facnrp = fnumsp.Pos(facet);

             ElementTransformation & trafom = ma->GetTrafo (eim, lh);
             ElementTransformation & trafop = ma->GetTrafo (eip, lh);

	     size_t tid= TaskManager::GetThreadId();
	     //auto nr = trace->StartTask(tid, tfac_getfe,
             //                           PajeTrace::Task::ID_TIMER);


             const FiniteElement & felm = fes.GetFE (eim, lh);
             const FiniteElement & felp = fes.GetFE (eip, lh);

	     //trace->StopTask(tid, nr);

             Array<int> dnumsm(felm.GetNDof(), lh);
             Array<int> dnumsp(felp.GetNDof(), lh);
             fes.GetDofNrs (eim, dnumsm);
             fes.GetDofNrs (eip, dnumsp);

             auto vnumsm = ma->GetElVertices (elm);
             auto vnumsp = ma->GetElVertices (elp);

             FlatVector<>
               elxm(dnumsm.Size()*num_equ, lh),
               elym(dnumsm.Size()*num_equ, lh),
               elxp(dnumsp.Size()*num_equ, lh),
               elyp(dnumsp.Size()*num_equ, lh);

             x.GetIndirect(dnumsm, elxm);
             x.GetIndirect(dnumsp, elxp);


             elym = 0.0;
             elyp = 0.0;

             int maxorder = max2 (felm.Order(), felp.Order());

             auto eltypem = trafom.GetElementType();
             auto eltypep = trafop.GetElementType();
             auto etfacet = ElementTopology::GetFacetType (eltypem, facnrm);

             Facet2ElementTrafo transformm(eltypem, vnumsm);
             Facet2ElementTrafo transformp(eltypep, vnumsp);

             SIMD_IntegrationRule ir_facet(etfacet, 2*maxorder);

             auto & ir_volm = transformm(facnrm, ir_facet, lh);
             auto & mirm = trafom(ir_volm, lh);

             auto & ir_volp = transformp(facnrp, ir_facet, lh);
             auto & mirp = trafop(ir_volp, lh);

             mirm.ComputeNormalsAndMeasure(eltypem, facnrm);

             // evaluate proxy-values
             ProxyUserData ud(numflux_trial_proxies.Size(), lh);
             const_cast<ElementTransformation&>(trafom).userdata = &ud;
             ud.fel = &felm;   // necessary to check remember-map
             // ud.elx = &elxm;
             // ud.lh = &lh;
             for (ProxyFunction * proxy : numflux_trial_proxies)
               ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);


	     //nr = trace->StartTask(tid, tfac_apply, PajeTrace::Task::ID_TIMER);

             for (ProxyFunction * proxy : numflux_trial_proxies)
               {
                 if (proxy->IsOther())
                   proxy->Evaluator()->Apply(felp, mirp, elxp,
                                             ud.GetAMemory(proxy));
                 else
                   proxy->Evaluator()->Apply(felm, mirm, elxm,
                                             ud.GetAMemory(proxy));
               }
	     //trace->StopTask(tid, nr);

             FlatMatrix<SIMD<double>> flux_pnt_values(num_equ, ir_facet.GetNIP(),
                                                      lh);
             numflux.Evaluate (mirm, flux_pnt_values);

             for (unsigned long i = 0; i < flux_pnt_values.Height(); i++)
               {
                 auto row = flux_pnt_values.Row(i);
                 for (size_t j = 0; j < row.Size(); j++)
                   row(j) *= mirm[j].GetWeight();
               }

	     //nr = trace->StartTask(tid, tfac_applyt, PajeTrace::Task::ID_TIMER);

	     evaluator->AddTrans(felm, mirm, flux_pnt_values, elym);
             evaluator->AddTrans(felp, mirp, flux_pnt_values, elyp);

	     // trace->StopTask(tid, nr);

             elym *= -1;
             y.AddIndirect(dnumsm, elym);
             y.AddIndirect(dnumsp, elyp);
           }
       });
  tfac.Stop();
}




class GradPhiCoefficientFunction : public CoefficientFunction
{
public:
  GradPhiCoefficientFunction (int _dim) : CoefficientFunction(_dim,false) { ; }


  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    throw Exception ("evaluate does nothing useful, needed since pure virtual");
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
			 BareSliceMatrix<SIMD<double>> values) const
  {
    MTP_UserData & ud = *(MTP_UserData*)ir.GetTransformation().userdata;

    for (size_t j = 0; j < ud.gradphi.Size(); j++)
      values.Row(j).AddSize(ir.Size()) = ud.gradphi(j);
  }
};





void ExportSymbolicDG ( py::module & m )
{
  m.def("ApplyDG",
        [] (shared_ptr<FESpace> fes,
            py::object Flux,
            py::object NumFlux,
            shared_ptr<BaseVector> x,
            shared_ptr<BaseVector> y)
        {

	  static Timer tcoef("ApplyDGPywrap");
	  RegionTimer reg(tcoef);

          py::object u = py::cast (fes).attr("TrialFunction")();
          py::object flux_u = Flux( u );
          // cout << "flux_u = " << flux_u << endl;
          py::object numflux_u = NumFlux( u, u.attr("Other")() );
          // cout << "numflux_u = " << numflux_u << endl;

          shared_ptr<CoefficientFunction> cpp_flux_u =
            py::extract<shared_ptr<CoefficientFunction>> (flux_u)();

          shared_ptr<CoefficientFunction> cpp_numflux_u =
            py::extract<shared_ptr<CoefficientFunction>> (numflux_u)();

	  cpp_numflux_u = Compile(cpp_numflux_u);
	  cpp_flux_u = Compile(cpp_flux_u);

          bool done = false;
          while (!done)
            {
              try
                {
                  ApplyDGOperator (*fes, *cpp_flux_u, *cpp_numflux_u,
                                   *x, *y);
                  done = true;
                }
              catch (LocalHeapOverflow &e)
                {
                  cout << "caught heap exception" << endl;
                  lh.CleanUp();
                  size_t newsize = 10*lh.Available();
                  // cout << "newsize = " << newsize << endl;
                  lh = LocalHeap (newsize, "applyDG", false);
                }
            }
        });

  m.def("GetTentGradPhi", [] (int dim) -> shared_ptr<CoefficientFunction>
        {
          return make_shared<GradPhiCoefficientFunction>(dim);
        });
}




