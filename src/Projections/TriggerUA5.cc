// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {


  void TriggerUA5::project(const Event& evt) {
    // Start with the assumption that the trigger fails
    _decision_sd = false;
    _decision_nsd = false;

    // Different trigger implementations for ppbar and pp!
    const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(evt, "CFS");
    foreach (const Particle& p, cfs.particles()) {
      const double eta = p.momentum().pseudorapidity();
      if (inRange(eta, -5.6, -2.0)) _n_minus++;
      else if (inRange(eta, 2.0, 5.6)) _n_plus++;
    }
    getLog() << Log::DEBUG << "Trigger -: " << _n_minus << ", Trigger +: " << _n_plus << endl;

    // Common SD/NSD trigger requirement: must activate at least one hodoscope
    if (_n_minus == 0 && _n_plus == 0) return; 
    _decision_sd = true;

    // Extra NSD trigger requirements
    const Beam& b = applyProjection<Beam>(evt, "Beam");
    _samebeams = (b.beams().first.pdgId() == b.beams().second.pdgId());
    if (_samebeams) {
      // PP
      if (_n_minus == 0 || _n_plus == 0) return;
    } else {
      // PPbar
      /// @todo Is this actually the exact trigger requirement?
      if (_n_minus * _n_plus < 4) return;
    }

    // Trigger success:
    _decision_nsd = true;
  }


}
