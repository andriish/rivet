// -*- C++ -*-
#include "Rivet/Analyses/CDF_2007_S7057202.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {

  
  const double CDF_2007_S7057202::_ybins[] = { 0.0, 0.1, 0.7, 1.1, 1.6, 2.1 };

  // Book histos and set counters for number of events passed in each one
  void CDF_2007_S7057202::init() {
    _histoD05 = bookHistogram1D(6, 1, 1, "0.1 < eta < 0.7, D=0.5");
    _histoD10 = bookHistogram1D(7, 1, 1, "0.1 < eta < 0.7, D=1.0");
    _binnedHistosD07.addHistogram(  0, 0.1, bookHistogram1D(1, 1, 1, "eta < 0.1, D=0.7"));
    _binnedHistosD07.addHistogram(0.1, 0.7, bookHistogram1D(2, 1, 1, "0.1 < eta < 0.7, D=0.7"));
    _binnedHistosD07.addHistogram(0.7, 1.1, bookHistogram1D(3, 1, 1, "0.7 < eta < 1.1, D=0.7"));
    _binnedHistosD07.addHistogram(1.1, 1.6, bookHistogram1D(4, 1, 1, "1.1 < eta < 1.6, D=0.7"));
    _binnedHistosD07.addHistogram(1.6, 2.1, bookHistogram1D(5, 1, 1, "1.6 < eta < 2.1, D=0.7"));

    size_t yind = 0;
    for (vector<AIDA::IHistogram1D*>::const_iterator histIt = _binnedHistosD07.getHistograms().begin();
        histIt != _binnedHistosD07.getHistograms().end(); ++histIt){
      _eventsPassed[*histIt] = 0.0;
      _yBinWidths[*histIt] = 2.0 * (_ybins[yind+1]-_ybins[yind]); 
      ++yind;
    }
    _eventsPassed[_histoD05] = 0.0;
    _yBinWidths[_histoD05] = 2.0*(-_ybins[1]+_ybins[2]);
    _eventsPassed[_histoD10] = 0.0;
    _yBinWidths[_histoD10] = 2.0*(-_ybins[1]+_ybins[2]);
  }


  // Do the analysis
  void CDF_2007_S7057202::analyze(const Event& event) {
    const double weight = event.weight();    
    
    const PseudoJets jetListD07 = applyProjection<FastJets>(event, "JetsD07").getPseudoJets();
    set< IHistogram1D*> passed;
    for (PseudoJets::const_iterator jet = jetListD07.begin(); jet != jetListD07.end(); ++jet) {
      const double pt = jet->perp();
      if (pt > _jetMinPT) {
        AIDA::IHistogram1D* histo = _binnedHistosD07.fill(fabs(jet->rapidity()), pt, weight);
        if (histo != 0) {
          if (histo->coordToIndex(pt) != IAxis::OVERFLOW_BIN) {
            passed.insert(histo);
	    _eventsPassed[histo] += weight;
          }
        }
      }
    }
    
    const PseudoJets jetListD05 = applyProjection<FastJets>(event, "JetsD05").getPseudoJets();
    for (PseudoJets::const_iterator jet = jetListD05.begin(); jet != jetListD05.end(); ++jet) {
      const double pt = jet->perp();
      if (pt > _jetMinPT) {
        double rap = fabs(jet->rapidity());
        if (rap >= _minY && rap < _maxY){
          _histoD05->fill(pt, weight);
          if (_histoD05->coordToIndex(pt) != IAxis::OVERFLOW_BIN){
            passed.insert(_histoD05);
	    _eventsPassed[_histoD05] += weight;
          }
        }
      }
    }
    
    const PseudoJets jetListD10 = applyProjection<FastJets>(event, "JetsD10").getPseudoJets();
    for (PseudoJets::const_iterator jet = jetListD10.begin(); jet != jetListD10.end(); ++jet){
      const double pt = jet->perp();
      if (pt > _jetMinPT) {
        double rap = fabs(jet->rapidity());
        if (rap >= _minY && rap < _maxY){
          _histoD10->fill(pt, weight);
          if (_histoD10->coordToIndex(pt) != IAxis::OVERFLOW_BIN){
            passed.insert(_histoD10);
	    _eventsPassed[_histoD10] += weight;
          }
        }
      }
    }
    
    //for (set<IHistogram1D*>::iterator histIt = passed.begin(); histIt != passed.end(); ++histIt) {
    // _eventsPassed[*histIt] += weight;
    //}
  }  



  // Normalise histograms to cross-section
  void CDF_2007_S7057202::finalize() {
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
