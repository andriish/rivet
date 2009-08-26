// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/E735_1998_S3905616.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  E735_1998_S3905616::E735_1998_S3905616()
      : Analysis("E735_1998_S3905616") {
    setBeams(PROTON, ANTIPROTON);
    const ChargedFinalState cfs;
    addProjection(cfs, "FS");
  }


  void E735_1998_S3905616::init() {
    _hist_multiplicity = bookHistogram1D(1, 1, 1);
  }


  void E735_1998_S3905616::analyze(const Event& event) {
    Log log = getLog();
    const ChargedFinalState& fs = applyProjection<ChargedFinalState>(event, "FS");
    const size_t numParticles = fs.particles().size();

    // Get the event weight
    const double weight = event.weight();

    // Fill histo of charged multiplicity distribution
    _hist_multiplicity->fill(numParticles, weight);
  }


  void E735_1998_S3905616::finalize() {
    normalize(_hist_multiplicity);
  }


}
