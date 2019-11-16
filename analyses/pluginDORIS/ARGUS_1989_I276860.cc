// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi+- K+- K0S pbar spectra continuum and Upsilon 1s 
  class ARGUS_1989_I276860 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1989_I276860);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // Book histograms
      _h_pi_ups1_p = bookHisto1D( 5, 1, 1);
      _h_pi_cont_p = bookHisto1D( 5, 1, 2);
      _h_pi_ups1_z = bookHisto1D( 9, 1, 1);
      _h_pi_cont_z = bookHisto1D( 9, 1, 2);
      
      _h_Kp_ups1_p = bookHisto1D( 6, 1, 1);
      _h_Kp_cont_p = bookHisto1D( 6, 1, 2);
      _h_Kp_ups1_z = bookHisto1D(10, 1, 1);
      _h_Kp_cont_z = bookHisto1D(10, 1, 2);
      
      _h_KS_ups1_p = bookHisto1D( 7, 1, 1);
      _h_KS_cont_p = bookHisto1D( 7, 1, 2);
      _h_KS_ups1_z = bookHisto1D(11, 1, 1);
      _h_KS_cont_z = bookHisto1D(11, 1, 2);
      
      _h_pt_ups1_p = bookHisto1D( 8, 1, 1);
      _h_pt_cont_p = bookHisto1D( 8, 1, 2);
      _h_pt_ups1_z = bookHisto1D(12, 1, 1);
      _h_pt_cont_z = bookHisto1D(12, 1, 2);

      // Counters
      _n_PiA[1] = bookCounter("/TMP/PiACont");
      _n_PiA[0] = bookCounter("/TMP/PiAUps1");
      _n_PiB[1] = bookCounter("/TMP/PiBCont");
      _n_PiB[0] = bookCounter("/TMP/PiBUps1");
      _n_Kp[1]  = bookCounter("/TMP/KpCont");
      _n_Kp[0]  = bookCounter("/TMP/KpUps1");
      _n_KS[1]  = bookCounter("/TMP/KSCont");
      _n_KS[0]  = bookCounter("/TMP/KSUps1");
      _n_ptA[1] = bookCounter("/TMP/ptACont");
      _n_ptA[0] = bookCounter("/TMP/ptAUps1");
      _n_ptB[1] = bookCounter("/TMP/ptBCont");
      _n_ptB[0] = bookCounter("/TMP/ptBUps1");
      // weights
      _weightSum_cont = 0.;
      _weightSum_Ups1 = 0.;
    }
    
   /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pdgId());
	if(id == 211 || id == 2212 || id == 310|| id==130 || id == 321) {
	  unstable.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      const FinalState & fs = apply<FinalState>(event, "FS");
      Particles upsilons = ufs.particles(Cuts::pid==553);
      const double weight = event.weight();
      // Continuum
      if (upsilons.empty()) { 
        MSG_DEBUG("No Upsilons found => continuum event");
        _weightSum_cont += weight;
	// stable particles
        foreach (const Particle& p, fs.particles(Cuts::abspid==211 or Cuts::abspid==2212 or Cuts::abspid==321 )) {
          const int id = abs(p.pdgId());
          const double xE = 2.*p.E()/sqrtS();
	  const double modp = p.p3().mod();
          const double beta = modp / p.E();
	  int idMom = !p.parents().empty() ? abs(p.parents()[0].pdgId()) : 0;
	  if(id==211) {
	    // not from K0S or Lambda decays
	    if(idMom!=310 && idMom != 3122) {
	      _h_pi_cont_p->fill(modp,weight);
	      _h_pi_cont_z->fill(xE  ,weight/beta);
	      _n_PiA[1]->fill(weight);
	    }
	    _n_PiB[1]->fill(weight);
	  }
	  else if(id==321) {
	    _h_Kp_cont_p->fill(modp,weight);
	    _h_Kp_cont_z->fill(xE  ,weight/beta);
	    _n_Kp[1]->fill(weight);
	  }
	  else if(id==2212) {
	    // not from K0S or Lambda decays
	    if(idMom!=310 && idMom != 3122) {
	      _h_pt_cont_p->fill(modp,weight);
	      _h_pt_cont_z->fill(xE  ,weight/beta);
	      _n_ptA[1]->fill(weight);
	    }
	    _n_ptB[1]->fill(weight);
	  }
	}
	// Unstable particles
        foreach (const Particle& p, ufs.particles(Cuts::pid==310 || Cuts::pid==130)) {
          const double xE = 2.*p.E()/sqrtS();
	  const double modp = p.p3().mod();
          const double beta = modp / p.E();
	  _h_KS_cont_p->fill(modp,weight);
	  _h_KS_cont_z->fill(xE  ,weight/beta);
	  _n_KS[1]->fill(weight);
	}
      }
      // Upsilon(s) found
      else { 
        MSG_DEBUG("Upsilons found => resonance event");
        for (const Particle& ups : upsilons) {
	  _weightSum_Ups1 += weight;
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups, unstable);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();
	  // loop over decay products
          foreach(const Particle& p, unstable) {
	    const int id = abs(p.pdgId());
            const FourMomentum p2 = cms_boost.transform(p.momentum());
	    const double xE = 2.*p2.E()/mass;
	    const double modp = p2.p3().mod();
	    const double beta = modp / p2.E();
	    int idMom = !p.parents().empty() ? abs(p.parents()[0].pdgId()) : 0;
	    if(id==211) {
	      // not from K0S or Lambda decays
	      if(idMom!=310 && idMom != 3122) {
		_h_pi_ups1_p->fill(modp,weight);
		_h_pi_ups1_z->fill(xE  ,weight/beta);
		_n_PiA[0]->fill(weight);
	      }
	      _n_PiB[0]->fill(weight);
	    }
	    else if(id==321) {
	      _h_Kp_ups1_p->fill(modp,weight);
	      _h_Kp_ups1_z->fill(xE  ,weight/beta);
	      _n_Kp[0]->fill(weight);
	    }
	    else if(id==2212) {
	      // not from K0S or Lambda decays
	      if(idMom!=310 && idMom != 3122) {
		_h_pt_ups1_p->fill(modp,weight);
		_h_pt_ups1_z->fill(xE  ,weight/beta);
		_n_ptA[0]->fill(weight);
	      }
	      _n_ptB[0]->fill(weight);
	    }
	    else if(id==310 || id==130) {
	      _h_KS_ups1_p->fill(modp,weight);
	      _h_KS_ups1_z->fill(xE  ,weight/beta);
	      _n_KS[0]->fill(weight);
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Scale histos
      if (_weightSum_cont > 0.) {
	scale(_h_pi_cont_p, 1./_weightSum_cont);
	scale(_h_Kp_cont_p, 1./_weightSum_cont);
	scale(_h_KS_cont_p, 1./_weightSum_cont);
	scale(_h_pt_cont_p, 1./_weightSum_cont);
	scale(_h_pi_cont_z, 1./_weightSum_cont);
	scale(_h_Kp_cont_z, 1./_weightSum_cont);
	scale(_h_KS_cont_z, 1./_weightSum_cont);
	scale(_h_pt_cont_z, 1./_weightSum_cont);

      }
      if (_weightSum_Ups1 > 0.) {
	scale(_h_pi_ups1_p, 1./_weightSum_Ups1);
	scale(_h_Kp_ups1_p, 1./_weightSum_Ups1);
	scale(_h_KS_ups1_p, 1./_weightSum_Ups1);
	scale(_h_pt_ups1_p, 1./_weightSum_Ups1);
	scale(_h_pi_ups1_z, 1./_weightSum_Ups1);
	scale(_h_Kp_ups1_z, 1./_weightSum_Ups1);
	scale(_h_KS_ups1_z, 1./_weightSum_Ups1);
	scale(_h_pt_ups1_z, 1./_weightSum_Ups1);
      }
      // Counters
      Scatter2DPtr nPiA = bookScatter2D(1, 1, 1, true);
      Scatter2DPtr nPiB = bookScatter2D(1, 1, 2, true);
      Scatter2DPtr nKS  = bookScatter2D(2, 1, 1, true);
      Scatter2DPtr nKp  = bookScatter2D(3, 1, 1, true);
      Scatter2DPtr nptA = bookScatter2D(4, 1, 1, true);
      Scatter2DPtr nptB = bookScatter2D(4, 1, 2, true);
      vector<double> scales = {_weightSum_Ups1,_weightSum_cont};
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_n_PiA[ix],1./scales[ix]);
	nPiA->point(ix).setY(_n_PiA[ix]->val(),_n_PiA[ix]->err());
	scale(_n_PiB[ix],1./scales[ix]);
	nPiB->point(ix).setY(_n_PiB[ix]->val(),_n_PiB[ix]->err());
	scale(_n_Kp[ix] ,1./scales[ix]);
	nKp->point(ix).setY(_n_Kp[ix]->val(),_n_Kp[ix]->err());
	scale(_n_KS[ix] ,1./scales[ix]);
	nKS->point(ix).setY(_n_KS[ix]->val(),_n_KS[ix]->err());
	scale(_n_ptA[ix],1./scales[ix]);
	nptA->point(ix).setY(_n_ptA[ix]->val(),_n_ptA[ix]->err());
	scale(_n_ptB[ix],1./scales[ix]);
	nptB->point(ix).setY(_n_ptB[ix]->val(),_n_ptB[ix]->err());
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pi_ups1_p, _h_pi_cont_p, _h_pi_ups1_z, _h_pi_cont_z;
    Histo1DPtr _h_Kp_ups1_p, _h_Kp_cont_p, _h_Kp_ups1_z, _h_Kp_cont_z;
    Histo1DPtr _h_KS_ups1_p, _h_KS_cont_p, _h_KS_ups1_z, _h_KS_cont_z; 
    Histo1DPtr _h_pt_ups1_p, _h_pt_cont_p, _h_pt_ups1_z, _h_pt_cont_z;

      // Counters
    CounterPtr _n_PiA[2], _n_PiB[2], _n_Kp[2], _n_KS[2], _n_ptA[2],_n_ptB[2];
    double _weightSum_cont,_weightSum_Ups1;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1989_I276860);

}
