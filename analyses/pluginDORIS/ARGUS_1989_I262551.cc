// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief phi spectrum in continuum, and Upsilon 1s and 2s decays
  class ARGUS_1989_I262551 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1989_I262551);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      _h_cont = bookHisto1D(1, 1, 1);
      _h_ups1 = bookHisto1D(2, 1, 1);
      _h_ups2 = bookHisto1D(2, 1, 2);
      _n_Phi[0]   = bookCounter("/TMP/NUps1");
      _n_Phi[1]   = bookCounter("/TMP/NUps2");
      _weightSum_cont = 0.;
      _weightSum_Ups1 = 0.;
      _weightSum_Ups2 = 0.;
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& phis) {
      for(const Particle & p: mother.children()) {
        const int id = p.pdgId();
	if(id == 333) {
	  phis.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, phis);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 or Cuts::pid==100553);
      const double weight = event.weight();
      // Continuum
      if (upsilons.empty()) { 
        MSG_DEBUG("No Upsilons found => continuum event");
        _weightSum_cont += weight;
        foreach (const Particle& p, ufs.particles(Cuts::pid==333)) {
          const double xp = 2.*p.E()/sqrtS();
          const double beta = p.p3().mod() / p.E();
	  _h_cont->fill(xp,weight/beta);
	}
      }
      // Upsilon(s) found
      else { 
        MSG_DEBUG("Upsilons found => resonance event");
        for (const Particle& ups : upsilons) {
          const int parentId = ups.pid();
	  if(parentId==553) {
	    _weightSum_Ups1 += weight;
	  }
	  else {
	    _weightSum_Ups2 += weight;
	  }
          Particles phis;
          // Find the decay products we want
          findDecayProducts(ups, phis);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();
	  // loop over decay products
          foreach(const Particle& p, phis) {
            const FourMomentum p2 = cms_boost.transform(p.momentum());
            const double xp = 2.*p2.E()/mass;
            const double beta = p2.p3().mod()/p2.E();
	    if(parentId==553) {
	      _n_Phi[0]->fill(weight);
	      _h_ups1->fill(xp,weight/beta);
	    }
	    else {
	      _n_Phi[1]->fill(weight);
	      _h_ups2->fill(xp,weight/beta);
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if (_weightSum_cont > 0.)
	scale(_h_cont, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      if (_weightSum_Ups1 > 0.) {
	scale(_h_ups1, 1./_weightSum_Ups1);
      }
      if (_weightSum_Ups2 > 0.) {
	scale(_h_ups2, 1./_weightSum_Ups2);
      }
      // Counters
      vector<double> scales = {_weightSum_Ups1,_weightSum_Ups2};
      for(unsigned int ix=0;ix<2;++ix) {
	Scatter2DPtr scatter = bookScatter2D(3+ix, 1, 1, true);
	scale(_n_Phi[ix],1./scales[ix]);
        scatter->point(0).setY(_n_Phi[ix]->val(),
			       _n_Phi[ix]->err());
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_cont, _h_ups1, _h_ups2;
    CounterPtr _n_Phi[2];
    double _weightSum_cont,_weightSum_Ups1,_weightSum_Ups2;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1989_I262551);


}
