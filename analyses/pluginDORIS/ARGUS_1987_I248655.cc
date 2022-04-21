// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief direct photon spectrum in upslion decay
  class ARGUS_1987_I248655 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1987_I248655);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      book(_h_gamma,2,1,1);
      book(_n_Ups,"TMP/nUps");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      for (const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==553)) {
	unsigned int nhadron(0);
	Particles photons;
	for(const Particle & child : ups.children()) {
	  if(PID::isHadron(child.pid()))
	    ++nhadron;
	  else if(child.pid()==PID::PHOTON)
	    photons.push_back(child);
	}
	if(nhadron!=0 && photons.empty()) {
	  _n_Ups->fill();
	}
	else if(nhadron!=0) {
	  LorentzTransform boost;
	  if (ups.p3().mod() > 1*MeV)
	    boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	  for(const Particle & gamma:photons)
	    _h_gamma->fill(2.*boost.transform(gamma.momentum()).E()/ups.mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_gamma, 1./ *_n_Ups);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_gamma;
    CounterPtr _n_Ups;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1987_I248655);

}
