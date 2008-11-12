// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2002_S4796047.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {


  // Book histograms
  void CDF_2002_S4796047::init() {
    _hist_multiplicity_630 = bookHistogram1D(1, 1, 1, "charged multiplicity at $\\sqrt{s}=630~\\text{GeV}$, $|\\eta| < 1$, $p_T>0.4~\\text{GeV}$");
    _hist_multiplicity_1800 = bookHistogram1D(2, 1, 1, "charged multiplicity at $\\sqrt{s}=1800~\\text{GeV}$, $|\\eta| < 1$, $p_T>0.4~\\text{GeV}$");
  }


  // Do the analysis
  void CDF_2002_S4796047::analyze(const Event& e) {
    Log log = getLog();

    const ChargedFinalState& fs = applyProjection<ChargedFinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Get the event weight
    const double weight = e.weight();

    // Get beams and average beam momentum
    const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
    const double meanBeamMom = ( beams.first.momentum().vector3().mod() + 
                                 beams.second.momentum().vector3().mod() ) / 2.0;

    if (2*meanBeamMom > 625 && 2*meanBeamMom < 635)
      _hist_multiplicity_630->fill(numParticles, weight);
    else if (2*meanBeamMom > 1795 && 2*meanBeamMom < 1805)
      _hist_multiplicity_1800->fill(numParticles, weight);

  }


  void CDF_2002_S4796047::finalize() {
    normalize(_hist_multiplicity_630, 3.21167);
    normalize(_hist_multiplicity_1800, 4.19121);
  }

}
