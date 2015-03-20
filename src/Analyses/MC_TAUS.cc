// -*- C++ -*-
#include "Rivet/Analyses/MC_ParticleAnalysis.hh"
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {


  /// @brief MC validation analysis for taus
  class MC_TAUS : public MC_ParticleAnalysis {
  public:

    MC_TAUS()
      : MC_ParticleAnalysis("MC_TAUS", 2, "tau")
    {    }


  public:

    void init() {
      IdentifiedFinalState taus;
      taus.acceptIdPair(PID::TAU);
      addProjection(taus, "Taus");

      MC_ParticleAnalysis::init();
    }


    void analyze(const Event& event) {
      const Particles taus = applyProjection<TauFinder>(event, "Taus").particlesByPt(0.5*GeV);
      MC_ParticleAnalysis::_analyze(event, taus);
    }


    void finalize() {
      MC_ParticleAnalysis::finalize();
    }

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TAUS);

}
