// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief OPAL photon/light meson paper
  /// @author Peter Richardson
  class OPAL_1998_S3749908 : public Analysis {
  public:

    /// Constructor
    OPAL_1998_S3749908()
      : Analysis("OPAL_1998_S3749908")
    {}


    /// @name Analysis methods
    //@{

    void init() {
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "FS");
      addProjection(UnstableFinalState(), "UFS");
      _histXePhoton   = bookHistogram1D( 2, 1, 1);
      _histXiPhoton   = bookHistogram1D( 3, 1, 1);
      _histXePi       = bookHistogram1D( 4, 1, 1);
      _histXiPi       = bookHistogram1D( 5, 1, 1);
      _histXeEta      = bookHistogram1D( 6, 1, 1);
      _histXiEta      = bookHistogram1D( 7, 1, 1);
      _histXeRho      = bookHistogram1D( 8, 1, 1);
      _histXiRho      = bookHistogram1D( 9, 1, 1);
      _histXeOmega    = bookHistogram1D(10, 1, 1);
      _histXiOmega    = bookHistogram1D(11, 1, 1);
      _histXeEtaPrime = bookHistogram1D(12, 1, 1);
      _histXiEtaPrime = bookHistogram1D(13, 1, 1);
      _histXeA0       = bookHistogram1D(14, 1, 1);
      _histXiA0       = bookHistogram1D(15, 1, 1);
    }


    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = applyProjection<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get event weight for histo filling
      const double weight = e.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.momentum().vector3().mod() +
                                   beams.second.momentum().vector3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");

      foreach (const Particle& p, ufs.particles()) {
        const int id = abs(p.pdgId());
	double xi = -log(p.momentum().vector3().mod()/meanBeamMom);
	double xE = p.momentum().t()/meanBeamMom;
        switch (id) {
	case 22:
	  _histXePhoton->fill(xE, weight);
	  _histXiPhoton->fill(xi, weight);
	  break;
	case 111:
	  _histXePi->fill(xE, weight);
	  _histXiPi->fill(xi, weight);
	  break;
	case 211:
	  _histXeEta->fill(xE, weight);
	  _histXiEta->fill(xi, weight);
	  break;
	case 213:
	  _histXeRho->fill(xE, weight);
	  _histXiRho->fill(xi, weight);
	  break;
	case 223:
	  _histXeOmega->fill(xE, weight);
	  _histXiOmega->fill(xi, weight);
	  break;
	case 331:
	  _histXeEtaPrime->fill(xE, weight);
	  _histXiEtaPrime->fill(xi, weight);
	  break;
	case 9000111:
	  _histXeA0->fill(xE, weight);
	  _histXiA0->fill(xi, weight);
	  break;
        }
      }
    }


    /// Finalize
    void finalize() {
      scale(_histXePhoton  , 1./sumOfWeights());
      scale(_histXiPhoton  , 1./sumOfWeights());
      scale(_histXePi      , 1./sumOfWeights());
      scale(_histXiPi      , 1./sumOfWeights());
      scale(_histXeEta     , 1./sumOfWeights());
      scale(_histXiEta     , 1./sumOfWeights());
      scale(_histXeRho     , 1./sumOfWeights());
      scale(_histXiRho     , 1./sumOfWeights());
      scale(_histXeOmega   , 1./sumOfWeights());
      scale(_histXiOmega   , 1./sumOfWeights());
      scale(_histXeEtaPrime, 1./sumOfWeights());
      scale(_histXiEtaPrime, 1./sumOfWeights());
      scale(_histXeA0      , 1./sumOfWeights());
      scale(_histXiA0      , 1./sumOfWeights());
    }

    //@}


  private:

      AIDA::IHistogram1D *_histXePhoton  ;
      AIDA::IHistogram1D *_histXiPhoton  ;
      AIDA::IHistogram1D *_histXePi      ;
      AIDA::IHistogram1D *_histXiPi      ;
      AIDA::IHistogram1D *_histXeEta     ;
      AIDA::IHistogram1D *_histXiEta     ;
      AIDA::IHistogram1D *_histXeRho     ;
      AIDA::IHistogram1D *_histXiRho     ;
      AIDA::IHistogram1D *_histXeOmega   ;
      AIDA::IHistogram1D *_histXiOmega   ;
      AIDA::IHistogram1D *_histXeEtaPrime;
      AIDA::IHistogram1D *_histXiEtaPrime;
      AIDA::IHistogram1D *_histXeA0      ;
      AIDA::IHistogram1D *_histXiA0      ;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1998_S3749908);

}
