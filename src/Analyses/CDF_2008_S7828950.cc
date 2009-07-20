// -*- C++ -*-
// CDF Run II inclusive jet cross-section using the midpoint algorithm.

#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S7828950.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  CDF_2008_S7828950::CDF_2008_S7828950()
    : Analysis("CDF_2008_S7828950")
  {
    setBeams(PROTON, ANTIPROTON);
    //setSqrtS(1960*GeV);
    const FinalState fs;
    addProjection(FastJets(fs, FastJets::CDFMIDPOINT, 0.7, 62.0*GeV), "JetsM07");
    setNeedsCrossSection(true);
  }

  
  const double CDF_2008_S7828950::_ybins[] = { 0.0, 0.1, 0.7, 1.1, 1.6, 2.1 };

  // Book histos and set counters for number of events passed in each one
  void CDF_2008_S7828950::init() {
    /// @todo What actually are these histos showing?
    _binnedHistosR07.addHistogram(  0, 0.1, bookHistogram1D(1, 1, 1));
    _binnedHistosR07.addHistogram(0.1, 0.7, bookHistogram1D(2, 1, 1));
    _binnedHistosR07.addHistogram(0.7, 1.1, bookHistogram1D(3, 1, 1));
    _binnedHistosR07.addHistogram(1.1, 1.6, bookHistogram1D(4, 1, 1));
    _binnedHistosR07.addHistogram(1.6, 2.1, bookHistogram1D(5, 1, 1));

    size_t yind = 0;
    foreach (AIDA::IHistogram1D* hist, _binnedHistosR07.getHistograms()) {
      _yBinWidths[hist] = 2.0 * (_ybins[yind+1]-_ybins[yind]);
      ++yind;
    }
  }


  // Do the analysis
  void CDF_2008_S7828950::analyze(const Event& event) {
    const double weight = event.weight();    
    
    foreach (const Jet& jet, applyProjection<FastJets>(event, "JetsM07").jets()) {
      _binnedHistosR07.fill(fabs(jet.momentum().rapidity()), jet.momentum().pT(), weight);
    }
  }  


  // Normalise histograms to cross-section
  void CDF_2008_S7828950::finalize() {
    foreach (AIDA::IHistogram1D* hist, _binnedHistosR07.getHistograms()) {
      scale(hist, crossSection()/nanobarn/sumOfWeights()/_yBinWidths[hist]);
    }
  }

  
}
