// -*- C++ -*-

// Field & Stuart underlying event analysis at CDF.
// Phys.Rev.D65:092002,2002 - no arXiv code.
// FNAL-PUB 01/211-E

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_MINBIAS.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {


  // Book histograms
  void CDF_2008_MINBIAS::init() {
    _hist_pt_vs_multiplicity = bookProfile1D(1, 1, 1, "Mean track pT vs multiplicity");
  }


  // Do the analysis
  void CDF_2008_MINBIAS::analyze(const Event& e) {
    Log log = getLog();

    const FinalState& fs = applyProjection<FinalState>(e, "FS");
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
//      const Vector3 mom3 = p->getMomentum().vector3();
//      const double energy = p->getMomentum().E();
      _hist_pt_vs_multiplicity->fill(numParticles, p->momentum().pT(), weight);
    }

  }


  void CDF_2008_MINBIAS::finalize() { 
  }


}
