// -*- C++ -*-
// CDF Run II inclusive jet cross-section using the midpoint algorithm.

#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S7828950.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  CDF_2008_S7828950::CDF_2008_S7828950()
    : _jetMinPT(62.0*GeV)
  {
    setBeams(PROTON, ANTIPROTON);
    //setSqrtS(1960*GeV);
    const FinalState fs;
    addProjection(FastJets(fs, FastJets::CDFMIDPOINT, 0.7), "JetsM07"); //??
    setNeedsCrossSection(true);
  }

  
  const double CDF_2008_S7828950::_ybins[] = { 0.0, 0.1, 0.7, 1.1, 1.6, 2.1 };

  // Book histos and set counters for number of events passed in each one
  void CDF_2008_S7828950::init() {
    /// @todo What actually are these histos showing?
    _binnedHistosR07.addHistogram(  0, 0.1, bookHistogram1D(1, 1, 1, "$\\eta < 0.1, R=0.7$",
							    "$E_{T}$ [GeV]","d^{2}\\sigma/dYdp_{T} [nb/GeV]"));
    _binnedHistosR07.addHistogram(0.1, 0.7, bookHistogram1D(2, 1, 1, "$0.1 < \\eta < 0.7, R=0.7$",
							    "$E_{T}$ [GeV]","d^{2}\\sigma/dYdp_{T} [nb/GeV]"));
    _binnedHistosR07.addHistogram(0.7, 1.1, bookHistogram1D(3, 1, 1, "$0.7 < \\eta < 1.1, R=0.7$",
							    "$E_{T}$ [GeV]","d^{2}\\sigma/dYdp_{T} [nb/GeV]"));
    _binnedHistosR07.addHistogram(1.1, 1.6, bookHistogram1D(4, 1, 1, "$1.1 < \\eta < 1.6, R=0.7$",
							    "$E_{T}$ [GeV]","d^{2}\\sigma/dYdp_{T} [nb/GeV]"));
    _binnedHistosR07.addHistogram(1.6, 2.1, bookHistogram1D(5, 1, 1, "$1.6 < \\eta < 2.1, R=0.7$",
							    "$E_{T}$ [GeV]","d^{2}\\sigma/dYdp_{T} [nb/GeV]"));

    size_t yind = 0;
    for (vector<AIDA::IHistogram1D*>::const_iterator histIt = _binnedHistosR07.getHistograms().begin();
        histIt != _binnedHistosR07.getHistograms().end(); ++histIt){
      _eventsPassed[*histIt] = 0.0;
      _yBinWidths[*histIt] = 2.0 * (_ybins[yind+1]-_ybins[yind]); 
      ++yind;
    }
  }


  // Do the analysis
  void CDF_2008_S7828950::analyze(const Event& event) {
    const double weight = event.weight();    
    
    const PseudoJets jetListM07 = applyProjection<FastJets>(event, "JetsM07").pseudoJets();
    set< IHistogram1D*> passed;
    for (PseudoJets::const_iterator jet = jetListM07.begin(); jet != jetListM07.end(); ++jet) {
      const double pt = jet->perp();
      if (pt > _jetMinPT) {
        AIDA::IHistogram1D* histo = _binnedHistosR07.fill(fabs(jet->rapidity()), pt, weight);
        if (histo != 0) {
          if (histo->coordToIndex(pt) != IAxis::OVERFLOW_BIN) {
            passed.insert(histo);
	    _eventsPassed[histo] += weight;
          }
        }
      }
    }    
  }  


  // Normalise histograms to cross-section
  void CDF_2008_S7828950::finalize() {
    Log log = getLog();
    const double xSecPerEvent = crossSection()/nanobarn / sumOfWeights();
    log << Log::INFO << "Cross-section = " << crossSection()/nanobarn << " nb" << endl;

    for (map<IHistogram1D*,double>::iterator histIt = _eventsPassed.begin(),
           histJt = _yBinWidths.begin(); histIt != _eventsPassed.end(); ++histIt, ++histJt) {
      IHistogram1D* hist = histIt->first;
      const double xSec = xSecPerEvent * histIt->second / histJt->second;
      normalize(hist, xSec);
    }
  }

  
}
