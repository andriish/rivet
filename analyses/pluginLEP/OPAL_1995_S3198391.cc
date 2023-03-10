// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief OPAL Delta++ fragmentation function paper
  ///
  /// @author Peter Richardson
  class OPAL_1995_S3198391 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_1995_S3198391);


    /// @name Analysis methods
    /// @{

    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_histXpDelta, 1, 1, 1);
    }


    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");

      for (const Particle& p : ufs.particles()) {
        if(p.abspid()==2224) {
          double xp = p.p3().mod()/meanBeamMom;
          _histXpDelta->fill(xp);
        }
      }
    }


    /// Finalize
    void finalize() {
      scale(_histXpDelta, 1./sumOfWeights());
    }

    /// @}


  private:

      Histo1DPtr _histXpDelta;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(OPAL_1995_S3198391, OPAL_1995_I398320);

}
