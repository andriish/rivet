// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief omega -> pi+pi-pi0
  class BESIII_2018_I1703033 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1703033);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_z  ,1,1,1);
      book(_h_phi,1,1,2);
      book(_h_cos,1,1,3);
      book(_h_ppi,1,1,4);
      book(_dalitz,"dalitz",50,-1,1,50,-1.1,0.9);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0) {
      for(const Particle & p : mother.children()) {
	int id = p.pid();
        if (id == PID::PIPLUS) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const std::complex<double> ii(0.,1.);
      for(const Particle& meson :
	    apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::OMEGA)) {
	// apply omega mass cut
	if(abs(meson.mass()-0.78265)>0.04) continue;
	// find decay products
      	unsigned int nstable(0);
       	Particles pip, pim, pi0;
      	findDecayProducts(meson, nstable, pip, pim, pi0);
	if(nstable !=3 || pi0.size()!=1 || pip.size()!=1 || pim.size()!=1) continue;
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(meson.momentum().betaVec());
	FourMomentum pp = boost.transform(pip[0].momentum());
	FourMomentum pm = boost.transform(pim[0].momentum());
	FourMomentum p0 = boost.transform(pi0[0].momentum());
	// variables from eqn 1 in paper (use version from Eqn 2.1 of 1010.3946)
	double Qc = meson.mass()-pi0[0].mass()-pim[0].mass()-pip[0].mass();
	double x = sqrt(3.)*(pm.E()-pp.E())/Qc;
	double y = 3./Qc*(p0.E()-pi0[0].mass())-1.;
	_dalitz->fill(x,y);
	// plot arg and phase of variable
	std::complex<double> zz = x+ii*y;
	_h_z->fill(norm(zz));
	double phi = arg(zz);
	if(phi<0.) phi+=2.*M_PI;
	_h_phi->fill(phi);
	double ctheta = -pp.p3().unit().dot(p0.p3().unit());
	_h_cos->fill(ctheta);
	_h_ppi->fill(p0.p3().mod());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_z  );
      normalize(_h_phi);
      normalize(_h_cos);
      normalize(_h_ppi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_z,_h_phi,_h_cos,_h_ppi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1703033);

}
