// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/DELPHI_1995_S3137023.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  // Constructor
  DELPHI_1995_S3137023::DELPHI_1995_S3137023() 
    : Analysis("DELPHI_1995_S3137023")
  {
    setBeams(ELECTRON, POSITRON); 
    addProjection(Beam(), "Beams");
    addProjection(ChargedFinalState(), "FS");
    addProjection(UnstableFinalState(), "UFS");
    _weightedTotalNumXiMinus = 0;
    _weightedTotalNumSigma1385Plus = 0;
  }


  void DELPHI_1995_S3137023::analyze(const Event& e) {
    // First, veto on leptonic events by requiring at least 4 charged FS particles
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
    if (numParticles < 2) {
      getLog() << Log::DEBUG << "Failed leptonic event cut" << endl;
      vetoEvent;
    }
    getLog() << Log::DEBUG << "Passed leptonic event cut" << endl;

    // Get event weight for histo filling
    const double weight = e.weight();

    // Get beams and average beam momentum
    const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
    const double meanBeamMom = ( beams.first.momentum().vector3().mod() + 
                                 beams.second.momentum().vector3().mod() ) / 2.0;
    getLog() << Log::DEBUG << "Avg beam momentum = " << meanBeamMom << endl;

    // Final state of unstable particles to get particle spectra
    const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");

    for (ParticleVector::const_iterator p = ufs.particles().begin(); p != ufs.particles().end(); ++p) {
      int id = abs(p->pdgId());
      switch (id) {
         case 3312:
            _histXpXiMinus->fill(p->momentum().vector3().mod()/meanBeamMom, weight);
            _weightedTotalNumXiMinus += weight;
            break;
         case 3114:
            _histXpSigma1385Plus->fill(p->momentum().vector3().mod()/meanBeamMom, weight);
            _weightedTotalNumSigma1385Plus += weight;
            break;
      }
    }

  }



  void DELPHI_1995_S3137023::init() {
    _histXpXiMinus       = bookHistogram1D(2, 1, 1);
    _histXpSigma1385Plus = bookHistogram1D(3, 1, 1);
  }

  // Finalize
  void DELPHI_1995_S3137023::finalize() { 
    normalize(_histXpXiMinus       , _weightedTotalNumXiMinus/sumOfWeights());
    normalize(_histXpSigma1385Plus , _weightedTotalNumSigma1385Plus/sumOfWeights());
  }



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<DELPHI_1995_S3137023> plugin_DELPHI_1995_S3137023;

}
