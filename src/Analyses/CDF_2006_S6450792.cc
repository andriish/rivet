// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2006_S6450792.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  CDF_2006_S6450792::CDF_2006_S6450792() : Analysis("CDF_2006_S6450792") {
    setBeams(PROTON, ANTIPROTON);
    setNeedsCrossSection(true);

    FinalState fs;
    addProjection(FastJets(fs, FastJets::CDFMIDPOINT, 0.7, 61.0*GeV), "ConeFinder");
  }


  void CDF_2006_S6450792::init() {

    _h_jet_pt = bookHistogram1D(1, 1, 1);

  }


  void CDF_2006_S6450792::analyze(const Event& event) {

    const Jets& jets = applyProjection<JetAlg>(event, "ConeFinder").jets();

    foreach (const Jet& jet, jets) {
      double y = fabs(jet.momentum().rapidity());
      if (y>0.1 && y<0.7) {
        _h_jet_pt->fill(jet.momentum().pT(), event.weight());
      }
    }
  
  }


  void CDF_2006_S6450792::finalize() {

    double delta_y = 1.2;
    scale(_h_jet_pt, crossSection()/nanobarn/sumOfWeights()/delta_y);

  }


}
