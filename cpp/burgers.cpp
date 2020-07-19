#include <solve.hpp>
using namespace ngsolve;

#include "conservationlaw.hpp"
#include "tconservationlaw_tp_impl.hpp"
#include <python_ngstd.hpp>

template <int D>
class Burgers : public T_ConservationLaw<Burgers<D>,D,1,1,false>
{
  Vec<D> b;
  typedef T_ConservationLaw<Burgers<D>,D,1,1,false> BASE;

public:
  Burgers (shared_ptr<TentPitchedSlab<D>> & atps, int order, const Flags & flags)
    : BASE (atps, order, flags) { ; }

  // these two were private
  using BASE::gfnu;
  using BASE::gfres;

  using BASE::CalcViscCoeffEl;
  using BASE::Flux;
  using BASE::u_reflect;
  using BASE::TransformBackIR;
  using BASE::CalcEntropy;


  // solve for û: Û = ĝ(x̂, t̂, û) - ∇̂ φ(x̂, t̂) ⋅ f̂(x̂, t̂, û)
  //
  // in this case, û = 2 Û / [ 1 + √(1 - 2 Û ⋅∇̂ φ(x̂, t̂)) ]
  //
  template <typename MIP=BaseMappedIntegrationPoint, typename SCAL=double>
  void TransformBack(const MIP & mip,
                     const Vec<D,SCAL > & grad,
                     const FlatVec<1,SCAL > u) const
  {
    SCAL sum = SCAL(0.0);
    for(size_t i : Range(D))
      sum += grad(i);
    u(0) = 2 * u(0) / (1.0 + sqrt(1.0 - 2.0*sum*u(0)) );
  }

  // solve for û at all points in an integration rule
  template <typename T>
  void TransformBackIR(const SIMD_BaseMappedIntegrationRule & mir,
                       FlatMatrix<T> grad, FlatMatrix<T> u) const
  {
    for (int i : Range(grad.Width()))
      {
	Vec<1,T> ui(u(0,i));
        Vec<D,T> gradi = grad.Col(i);
        TransformBack<SIMD<BaseMappedIntegrationPoint>, T>(mir[i],gradi,ui);
        u(0,i) = ui(0);
      }
  }

  // Flux on element
  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    for(size_t i : Range(mir))
      flux.Col(i) = 0.5 * u(0,i)*u(0,i);
  }

  // Helper for Flux on element (compute a column)
  template<typename SCAL=double>
  Vec<1,SCAL> Flux (const FlatVec<1,SCAL> & ul, const FlatVec<1,SCAL> & ur,
                    const Vec<D,SCAL> & nv) const
  {
    SCAL sumn = SCAL(0.0);
    for(size_t i : Range(D))
      sumn += nv(i);
    // solution from inside, in case of a boundary facet
    SCAL fln = ul(0)*ul(0)/2.0 * sumn;
    SCAL frn = ur(0)*ur(0)/2.0 * sumn;
    SCAL um = 0.5*(ul(0)+ur(0));

    SCAL fprimen = um * sumn; // derivative of f with respect to u

    SCAL fluxn = IfPos(fprimen, fln, frn);

    return Vec<1,SCAL>(fluxn);
  }

  // Numerical Flux on a facet.  ul and ur are the values at integration points
  // of the two elements adjacent to an internal facet
  // of the spatial mesh of a tent.
  void Flux(FlatMatrix<SIMD<double>> ula, FlatMatrix<SIMD<double>> ura,
            FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    for(size_t i : Range(ula.Width()))
      {
        SIMD<double> sumn = 0.0;
        for(size_t j : Range(D))
          sumn += normals(j,i);

        SIMD<double> um = 0.5*(ula(0,i) + ura(0,i));
        SIMD<double> fprimen = um * sumn; // f'(u) * n gives characteristic speed
        fna(0,i) = IfPos(fprimen,
                         0.5*(ula(0,i)*ula(0,i))*sumn,
                         0.5*(ura(0,i)*ura(0,i))*sumn);
      }
  }

  void u_reflect(FlatMatrix<SIMD<double>> u,
                 FlatMatrix<SIMD<double>> normals,
                 FlatMatrix<SIMD<double>> u_refl) const
  {
    cout << "no reflecting B.C. for Burgers Equation" << endl;
  }

  void CalcEntropy(FlatMatrix<AutoDiff<1,SIMD<double>>> adu,
                   FlatMatrix<AutoDiff<1,SIMD<double>>> grad,
		   FlatMatrix<SIMD<double>> dEdt, FlatMatrix<SIMD<double>> F) const
  {
    // E = u^2/2,
    // F = u^3/3 [ 1 ]*D
    for(size_t i : Range(adu.Width()))
      {
        auto ui = adu(0,i);
        auto sumgrad = grad(0,i);
        auto Fi = ui*ui*ui/3.0;
        for(size_t j : Range(1,D))
          sumgrad += grad(j,i);

        AutoDiff<1,SIMD<double>> adE = ui*ui/2.0 - Fi*sumgrad; // E - (F,grad(phi))
        dEdt(0,i) = adE.DValue(0);
        for(size_t j : Range(D))
          F(j,i) = Fi.Value();
      }
  }


  void EntropyFlux (FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
                    FlatMatrix<SIMD<double>> n,
                    FlatMatrix<SIMD<double>> flux) const
  {
    for(size_t i : Range(ul.Width()))
      {
        SIMD<double> sumn = 0.0;
        for(size_t j : Range(D))
          sumn += n(j,i);
	//////// upwind flux ////////
	auto um = 0.5 * (ul(0,i) + ur(0,i));
	auto Fprimen = um*um * sumn;
	flux(0,i) = IfPos(Fprimen,
			  1./3*ul(0,i)*ul(0,i)*ul(0,i)*sumn,
			  1./3*ur(0,i)*ur(0,i)*ur(0,i)*sumn);
      }
  }

  // Given the entropy, entropy residual and 'hi' in the cylinder
  // at points in a MIR, where
  // hi = ([Jacobian determinant for element]/DIM)^{1/DIM}
  //
  // compute the viscosity coefficient
  // ν = max_{T ∈  Tᵢ} min(ν_*ᵀ, νₑᵀ), where
  // ν_*ᵀ = κ₂ diam(T)||Dᵤ f̂||_{L^∞(T)},
  // vₑᵀ = c_X² ||Rₕ||_{L^∞ (T)} / |E̅|, where E̅ is the mean value of entropy
  // c_X² is an effective local grid size of X, Rₕ is entropy residual
  //
  // betai corresponds to v_*ᵀ and visci * hi/Emean corresponds to vₑᵀ.
  void CalcViscCoeffEl(const SIMD_BaseMappedIntegrationRule & mir,
                       FlatMatrix<SIMD<double>> elu_ipts,
                       FlatMatrix<SIMD<double>> res_ipts,
                       const double hi, double & coeff) const
  {
    int nipt = mir.IR().GetNIP();

    double betai = 0.0;
    double visci = 0.0;
    double Emean = 0.0;

    for(size_t i : Range(elu_ipts.Width()))
      {
        auto res = res_ipts(0,i);
        auto ui = elu_ipts(0,i);
        for(size_t j : Range(SIMD<double>::Size()))
          {
            if(fabs(res[j]) > visci)
              visci = fabs(res[j]);
            Emean += ui[j]*ui[j]/2.0;
            betai = max(betai, fabs(ui[j]));
          }
      }
    Emean /= nipt;
    coeff = min (visci * hi/Emean, betai);
    coeff *= hi;
  }
};

/////////////////////////////////////////////////////////////////////////
//                 EXPORT TO PYTHON

typedef CoefficientFunction CF;
typedef GridFunction GF;
typedef Burgers<1> B1;
typedef Burgers<2> B2;

void ExportConsLaw(py::module & m)
{
  py::class_<B1, shared_ptr<B1>>
    (m,"Burgers1", "Burgers equation in 1 spatial dimension")
    .def(py::init([](shared_ptr<TentPitchedSlab<1>> & tps,
                     int order, py::dict pyflags)
                  {
                    const Flags flags = py::extract<Flags> (pyflags)();
                    auto cl = make_shared<B1>(tps, order, flags);
                    cl->SetBC(); //use old style bc numbers for now
                    return cl;
                  }),
         py::arg("tentslab"),
         py::arg("order"),
         py::arg("flags")=py::dict())
    .def_property_readonly("sol", [](shared_ptr<B1> self) { return self->gfu; })
    .def_property_readonly("res", [](shared_ptr<B1> self) { return self->gfres; })
    .def_property_readonly("nu", [](shared_ptr<B1> self) { return self->gfnu; })
    .def_property_readonly("space",
                           [](shared_ptr<B1> self) -> shared_ptr<FESpace>
                           { return self->fes; })

    // Set the initial data
    .def("SetInitial",
         [](shared_ptr<B1> self, shared_ptr<CF> cf)
         {
           SetValues(cf,*(self->gfu),VOL,0,*(self->pylh));
           self->uinit = self->u; // set data used for b.c.
         })
    .def("PropagatePicard",
         [](shared_ptr<B1> self, shared_ptr<BaseVector> vecu, int steps)
         {
           if(steps==-1)
             steps = 1;
           self->PropagatePicard(steps,*vecu,self->uinit,*(self->pylh));
         }, py::arg("vec"),py::arg("steps")=-1)
    .def("Tau",[](shared_ptr<B1> self) { return self->gftau; })

    ; // please keep me on my own line

  py::class_<B2, shared_ptr<B2>>
    (m,"Burgers2", "Burgers equation in 2 spatial dimensions")
    .def(py::init([](shared_ptr<TentPitchedSlab<2>> & tps,
                     int order, py::dict pyflags)
                  {
                    const Flags flags = py::extract<Flags> (pyflags)();
                    auto cl = make_shared<B2>(tps, order, flags);
                    cl->SetBC(); //use old style bc numbers for now
                    return cl;
                  }),
         py::arg("tentslab"),
         py::arg("order"),
         py::arg("flags")=py::dict())
    .def_property_readonly("sol", [](shared_ptr<B2> self) { return self->gfu; })
    .def_property_readonly("res", [](shared_ptr<B2> self) { return self->gfres; })
    .def_property_readonly("nu", [](shared_ptr<B2> self) { return self->gfnu; })
    .def_property_readonly("space",
                           [](shared_ptr<B2> self) -> shared_ptr<FESpace>
                           { return self->fes; })

    // Set the initial data
    .def("SetInitial",
         [](shared_ptr<B2> self, shared_ptr<CF> cf)
         {
           SetValues(cf,*(self->gfu),VOL,0,*(self->pylh));
           self->uinit = self->u; // set data used for b.c.
         })
    .def("PropagatePicard",
         [](shared_ptr<B2> self, shared_ptr<BaseVector> vecu, int steps)
         {
           if(steps==-1)
             steps = 1;
           self->PropagatePicard(steps,*vecu,self->uinit,*(self->pylh));
         }, py::arg("vec"),py::arg("steps")=-1)
    .def("Tau",[](shared_ptr<B2> self) { return self->gftau; })

    ; // please keep me on my own line

}
