// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> X_s gamma spectrum
  class BABAR_2012_I1123662 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1123662);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_spectrum, 1, 1, 1);
      book(_nBottom, "TMP/BottomCounter");
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over bottoms
      for (const Particle& bottom : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==521 ||
										     Cuts::abspid==511)) {
	// remove mixing entries etc
	for (const Particle & child : bottom.children())
          if (child.abspid() == 511 || child.pid()==bottom.pid() ) continue;
        _nBottom->fill();
	FourMomentum pgamma(0.,0.,0.,0.);
	unsigned int ngamma = 0;
        for (const Particle & child : bottom.children()) {
	  if (child.pid() == PID::PHOTON) {
            ngamma += 1;
            pgamma += child.momentum();
          }
	}
	if (ngamma != 1) continue;
	const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(bottom.momentum().betaVec());
	double eGamma = boost.transform(pgamma).E();
	_h_spectrum->fill(eGamma);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spectrum, 1e5/_nBottom->sumW());
    }    

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spectrum;
    CounterPtr _nBottom;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1123662);

}
