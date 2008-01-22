// -*- C++ -*-
#include "Rivet/Analyses/CDF_2007_S7057202.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {
  
  // Book histos and set counters for number of events passed in each one
  void CDF_2007_S7057202::init() {
    
    _eventsTried = 0.0;

    /// @todo XML escaping should be done by LWH
    _histoD05 = bookHistogram1D(6,1,1,"0.1 &lt; eta &lt; 0.7, D=0.5");
    _histoD10 = bookHistogram1D(7,1,1,"0.1 &lt; eta &lt; 0.7, D=1.0");

    _binnedHistosD07.addHistogram(  0, 0.1, bookHistogram1D(1,1,1,"eta &lt; 0.1, D=0.7"));
    _binnedHistosD07.addHistogram(0.1, 0.7, bookHistogram1D(2,1,1,"0.1 &lt; eta &lt; 0.7, D=0.7"));
    _binnedHistosD07.addHistogram(0.7, 1.1, bookHistogram1D(3,1,1,"0.7 &lt; eta &lt; 1.1, D=0.7"));
    _binnedHistosD07.addHistogram(1.1, 1.6, bookHistogram1D(4,1,1,"1.1 &lt; eta &lt; 1.6, D=0.7"));
    _binnedHistosD07.addHistogram(1.6, 2.1, bookHistogram1D(5,1,1,"1.6 &lt; eta &lt; 2.1, D=0.7"));
  
    for(vector<AIDA::IHistogram1D*>::const_iterator histIt = _binnedHistosD07.getHistograms().begin();
        histIt != _binnedHistosD07.getHistograms().end(); ++histIt){
      _eventsPassed[*histIt] = 0.0;
    }
    
    _eventsPassed[_histoD05] = 0.0;
    _eventsPassed[_histoD10] = 0.0;

  }

////////////////////////////////////////////////////////////////////////////////
  void CDF_2007_S7057202::analyze(const Event& event) {
    const double weight = event.weight();    
    event.applyProjection(_ktprojD07);
    event.applyProjection(_ktprojD05);
    event.applyProjection(_ktprojD10);
    
    PseudoJets jetList = _ktprojD07.getPseudoJets();
    
    set< IHistogram1D*> passed;
    for (PseudoJets::const_iterator jet = jetList.begin(); 
         jet != jetList.end(); ++jet) {
      const double pt = jet->perp();
      if (pt > _jetMinPT) {
        
        AIDA::IHistogram1D* histo = 
        _binnedHistosD07.fill(fabs(jet->rapidity()), pt, weight);
        
        if(histo != 0){
          if(histo->coordToIndex(pt) != IAxis::OVERFLOW_BIN){
            passed.insert(histo);
          }
        }
      }
    }
    
    PseudoJets jetListD05 = _ktprojD05.getPseudoJets();
    for (PseudoJets::iterator jet = jetListD05.begin(); 
         jet != jetListD05.end(); ++jet) {
      const double pt = jet->perp();
      if (pt > _jetMinPT) {
        double rap = fabs(jet->rapidity());
        if (rap >= _minY && rap < _maxY){
          _histoD05->fill(pt, weight);
          if(_histoD05->coordToIndex(pt) != IAxis::OVERFLOW_BIN){
            passed.insert(_histoD05);
          }
        }
      }
    }
    
    PseudoJets jetListD10 = _ktprojD10.getPseudoJets();
    for(PseudoJets::iterator jet = jetListD10.begin(); 
        jet != jetListD10.end(); ++jet){
      const double pt = jet->perp();
      if (pt > _jetMinPT) {
        double rap = fabs(jet->rapidity());
        if (rap >= _minY && rap < _maxY){
          _histoD10->fill(pt, weight);
          if(_histoD10->coordToIndex(pt) != IAxis::OVERFLOW_BIN){
            passed.insert(_histoD10);
          }
        }
      }
    }
    
    // Increment the event counters for each histogram
    _eventsTried += weight;
    
    for (set<IHistogram1D*>::iterator histIt = passed.begin(); 
         histIt != passed.end(); ++histIt) {      
      _eventsPassed[*histIt] += weight;
    }
    return;
  }  
//////////////////////////////////////////////////////////////////////////////// 
  // Normalise histograms to cross-section
  void CDF_2007_S7057202::finalize() {
    double xSecPerEvent = crossSection() / _eventsTried;
    // HepData data is in nb, crossSection returns pb.
    /// @todo Choose consistent units set
    xSecPerEvent = 0.001 * xSecPerEvent; 
    
    for (map<IHistogram1D*,double>::iterator histIt = _eventsPassed.begin();
         histIt != _eventsPassed.end(); ++histIt) {
      IHistogram1D* hist = histIt->first;
      double xSec = xSecPerEvent * histIt->second;
      normalize(hist, xSec);
    }
    return;
  }
  
}
////////////////////////////////////////////////////////////////////////////////