// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief k0s in b decays
  class ARGUS_1994_I354224 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1994_I354224);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _h_K = bookHisto1D(1, 1, 1);
      _nB=0.;
    }

    void analyzeDecay(Particle mother, Particles & kaons) {
      for(const Particle & p : mother.children()) {
	if(p.pdgId()==310) {
	  kaons.push_back(p);
	}
	else if(!p.children().empty()) {
	  analyzeDecay(p,kaons);
	} 
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==511 or Cuts::abspid==521)) {
	if(!p.children().empty()) {
	  if(p.children()[0].pdgId()==p.pdgId()) continue;
	}
	FourMomentum pB = p.momentum();
	const LorentzTransform B_boost = LorentzTransform::mkFrameTransformFromBeta(pB.betaVec());
	_nB += event.weight();
	Particles kaons;
	analyzeDecay(p,kaons);
	for(const Particle & kaon : kaons) {
	  FourMomentum pKaon = B_boost.transform(kaon.momentum());
	  _h_K->fill(pKaon.p3().mod(),event.weight());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_K, 1./_nB);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_K;
    double _nB;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1994_I354224);


}
