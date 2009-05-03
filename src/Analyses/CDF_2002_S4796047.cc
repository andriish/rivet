// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2002_S4796047.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  CDF_2002_S4796047::CDF_2002_S4796047()
  { 
    setBeams(PROTON, ANTIPROTON);
    addProjection(Beam(), "Beams");
    const ChargedFinalState cfs(-1.0, 1.0, 0.4*GeV);
    addProjection(cfs, "FS");
  }


  // Book histograms
  void CDF_2002_S4796047::init() {
    /// @todo Cross-section units
    _hist_multiplicity_630 = bookHistogram1D(1, 1, 1, 
                                             "Charged multiplicity at $\\sqrt{s} = 630~\\text{GeV}$, $|\\eta| < 1$, $p_T > 0.4~\\text{GeV}$",
                                             "$n_\\text{ch}$", "\\d{\\sigma}/\\d{n_\\text{ch}}$");
    /// @todo Cross-section units
    _hist_multiplicity_1800 = bookHistogram1D(2, 1, 1, 
                                              "Charged multiplicity at $\\sqrt{s} = 1800~\\text{GeV}$, $|\\eta| < 1$, $p_T > 0.4~\\text{GeV}$",
                                              "$n_\\text{ch}$", "\\d{\\sigma}/\\d{n_\\text{ch}}$");
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
    const double sumBeamMom = ( beams.first.momentum().vector3().mod() + 
                                beams.second.momentum().vector3().mod() );

    if (inRange(sumBeamMom/GeV, 625, 635)) {
      _hist_multiplicity_630->fill(numParticles, weight);
    } else if (inRange(sumBeamMom/GeV, 1795, 1805)) {
      _hist_multiplicity_1800->fill(numParticles, weight);
    }
  }


  void CDF_2002_S4796047::finalize() {
    /// @todo Get cross-section from the generator
    normalize(_hist_multiplicity_630, 3.21167);
    normalize(_hist_multiplicity_1800, 4.19121);
  }


}
