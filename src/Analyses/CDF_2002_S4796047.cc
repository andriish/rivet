// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2002_S4796047.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  CDF_2002_S4796047::CDF_2002_S4796047()
    : Analysis("CDF_2002_S4796047")
  { 
    setBeams(PROTON, ANTIPROTON);
    addProjection(Beam(), "Beam");
    const ChargedFinalState cfs(-1.0, 1.0, 0.4*GeV);
    addProjection(cfs, "FS");
  }


  // Book histograms
  void CDF_2002_S4796047::init() {
    /// @todo Cross-section units
    _hist_multiplicity_630  = bookHistogram1D(1, 1, 1);
    /// @todo Cross-section units
    _hist_multiplicity_1800 = bookHistogram1D(2, 1, 1);
    _hist_pt_vs_multiplicity_630  = bookProfile1D(3, 1, 1);
    _hist_pt_vs_multiplicity_1800 = bookProfile1D(4, 1, 1);
  }


  // Do the analysis
  void CDF_2002_S4796047::analyze(const Event& e) {
    Log log = getLog();
    const double sqrtS = applyProjection<Beam>(e, "Beam").sqrtS();
    const ChargedFinalState& fs = applyProjection<ChargedFinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Get the event weight
    const double weight = e.weight();

    // Fill histos of charged multiplicity distributions
    if (fuzzyEquals(sqrtS, 630/GeV)) {
      _hist_multiplicity_630->fill(numParticles, weight);
    } 
    else if (fuzzyEquals(sqrtS, 1800/GeV)) {
      _hist_multiplicity_1800->fill(numParticles, weight);
    }

    // Fill histos for <pT> vs. charged multiplicity
    foreach (const Particle& p, fs.particles()) {
      const double pT = p.momentum().pT() / GeV;

      if (fuzzyEquals(sqrtS, 630/GeV)) {
        _hist_pt_vs_multiplicity_630->fill(numParticles, pT, weight);
      }
      else if (fuzzyEquals(sqrtS, 1800/GeV)) {
        _hist_pt_vs_multiplicity_1800->fill(numParticles, pT, weight);
      }
    }

  }


  void CDF_2002_S4796047::finalize() {
    /// @todo Get cross-section from the generator
    normalize(_hist_multiplicity_630, 3.21167);
    normalize(_hist_multiplicity_1800, 4.19121);
  }


}
