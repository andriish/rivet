// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_NOTE_9337.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {


  // Book histograms
  void CDF_2008_NOTE_9337::init() {
    _hist_pt_vs_multiplicity = bookProfile1D(1, 1, 1, "Mean track $p_\\perp$ vs multiplicity");
  }


  // Do the analysis
  void CDF_2008_NOTE_9337::analyze(const Event& e) {
    Log log = getLog();

    const ChargedFinalState& fs = applyProjection<ChargedFinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
    if (numParticles < 1) {
      getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
      vetoEvent(e);
    }

    // Get the event weight
    const double weight = e.weight();

    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      // Get momentum and energy of each particle.
//      const Vector3 mom3 = p->momentum().vector3();
//      const double energy = p->momentum().E();
      _hist_pt_vs_multiplicity->fill(numParticles, p->momentum().pT(), weight);
    }

  }


  void CDF_2008_NOTE_9337::finalize() { 
  }


}
