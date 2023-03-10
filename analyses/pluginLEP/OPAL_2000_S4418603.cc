// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief OPAL K0 fragmentation function paper
  ///
  /// @author Peter Richardson
  class OPAL_2000_S4418603 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_2000_S4418603);


    /// @name Analysis methods
    /// @{

    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_histXeK0   , 3, 1, 1);
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
        const int id = p.abspid();
        if (id == PID::K0S || id == PID::K0L) {
          double xE = p.E()/meanBeamMom;
          _histXeK0->fill(xE);
        }
      }
    }


    /// Finalize
    void finalize() {
      scale(_histXeK0, 1./sumOfWeights());
    }

    /// @}


  private:

    /// Histogram
    Histo1DPtr _histXeK0;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(OPAL_2000_S4418603, OPAL_2000_I529898);

}
