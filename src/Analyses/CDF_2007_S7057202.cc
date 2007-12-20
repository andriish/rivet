// -*- C++ -*-
#include "Rivet/Analyses/CDF_2007_S7057202.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {
  
  //const double CDF_2007_S7057202::_ybins[] = {0.1, 0.7, 1.1, 1.6, 2.1};

  
  // Book histos and set counters for number of events passed in each one
  void CDF_2007_S7057202::init() {
    _eventsTried = 0.0;

    /// @todo Indexing a map by double is a bad idea...
    _histosD07[_ybins[0]] = bookHistogram1D(1, 1, 1,"eta < 0.1, D = 0.7");
    _histosD07[_ybins[1]] = bookHistogram1D(2, 1, 1,"0.1 < eta < 0.7, D = 0.7");
    _histosD07[_ybins[2]] = bookHistogram1D(3, 1, 1,"0.7 < eta < 1.1, D = 0.7");
    _histosD07[_ybins[3]] = bookHistogram1D(4, 1, 1,"1.1 < eta < 1.6, D = 0.7");
    _histosD07[_ybins[4]] = bookHistogram1D(5, 1, 1,"1.6 < eta < 2.1, D = 0.7");
    _histoD05 = bookHistogram1D(6, 1, 1,"0.1 < eta < 0.7, D = 0.5");
    _histoD10 = bookHistogram1D(7, 1, 1,"0.1 < eta < 0.7, D = 1.0");
    
    for (map<double, IHistogram1D*>::iterator histIt = _histosD07.begin();
         histIt != _histosD07.end(); ++histIt) {
      _eventsPassedD07[histIt->second] = 0.0;
    }
    
    _eventsPassedD05 = 0.0;
    _eventsPassedD10 = 0.0;
  }

  
  void CDF_2007_S7057202::analyze(const Event& event) {
    const double weight = event.weight();    
    event.applyProjection(_ktprojD07);
    event.applyProjection(_ktprojD05);
    event.applyProjection(_ktprojD10);
    
    /// @todo This is a pretty fiddly way to fill the appropriate histo...
    PseudoJets jetList = _ktprojD07.getPseudoJets();
    set<IHistogram1D*> passed;
    for (PseudoJets::const_iterator jet = jetList.begin(); jet != jetList.end(); ++jet) {
      const double pt = jet->perp();
      
      if (pt > _jetMinPT) {
        map<double, IHistogram1D*>::iterator histIt = _histosD07.upper_bound(jet->eta());
        
        if(histIt != _histosD07.end()){
          IHistogram1D* histo = histIt->second;
          // ...fill the histo
          histo->fill(pt, weight);
          if(histo->coordToIndex(pt) != IAxis::OVERFLOW_BIN){
            passed.insert(histo);
          }
        }
      }
    }
    
    // Increment the event counters for each histogram
    _eventsTried += weight;
    
    for (set<IHistogram1D*>::iterator histIt = passed.begin(); histIt != passed.end(); ++histIt) {      
      _eventsPassedD07[*histIt] += weight;
    }
    

    // Do the same for the D=0.5 jets
    event.applyProjection(_ktprojD05);
    PseudoJets jetListD05 = _ktprojD07.getPseudoJets();
    for (PseudoJets::iterator jet = jetListD05.begin(); jet != jetListD05.end(); ++jet) {
      const double pt = jet->perp();
      if (pt > _jetMinPT) {
        double rap = fabs(jet->rapidity());
        if (rap >= _ybins[0] && rap < _ybins[1]) _histoD05->fill(pt, weight);
      }
    }
    _eventsPassedD05 += weight;
    
    // Do the same for the D=1.0 jets
    event.applyProjection(_ktprojD10);
    PseudoJets jetListD10 = _ktprojD10.getPseudoJets();
    for(PseudoJets::iterator jet = jetListD10.begin(); jet != jetListD10.end(); ++jet){
      const double pt = jet->perp();
      if (pt > _jetMinPT) {
        double rap = fabs(jet->rapidity());
        if (rap >= _ybins[0] && rap < _ybins[1]) _histoD10->fill(pt, weight);
      }
    }
    _eventsPassedD10 += weight;
    
  }
  


 
  // Normalise histograms to cross-section
  void CDF_2007_S7057202::finalize() {
    double xSecPerEvent = crossSection() / _eventsTried;
    // HepData data is in nb, crossSection returns pb.
    /// @todo Choose consistent units set
    xSecPerEvent = 0.001 * xSecPerEvent; 
    
    /// @todo Check diff between this and the ones below.
    for (map<IHistogram1D*,double>::iterator histIt = _eventsPassedD07.begin();
         histIt != _eventsPassedD07.end(); ++histIt) {
      IHistogram1D* hist = histIt->first;
      double xSec = xSecPerEvent * histIt->second;
      /// @todo Can't the normalize() function be used here?
      int nBins = hist->axis().bins();
      double hArea = 0.0;
      for(int iBin = 0; iBin != nBins; ++iBin){
        hArea += hist->binHeight(iBin) * hist->axis().binWidth(iBin);
      }
      hist->scale(xSec / hArea);
    }
    
    // Do the same for D05 histogram
    {
      const double xSec = xSecPerEvent * _eventsPassedD05;
      const size_t nBins = _histoD05->axis().bins();
      double hArea = 0.0;
      for (size_t iBin = 0; iBin != nBins; ++iBin) {
        hArea += _histoD05->binHeight(iBin) * _histoD05->axis().binWidth(iBin);
      }
      _histoD05->scale(xSec / hArea);
    }
    
    // Do the same for D10 histogram
    {
      const double xSec = xSecPerEvent * _eventsPassedD10;
      const size_t nBins = _histoD10->axis().bins();
      double hArea = 0.0;
      for (size_t iBin = 0; iBin != nBins; ++iBin) {
        hArea += _histoD10->binHeight(iBin) * _histoD10->axis().binWidth(iBin);
      }
      _histoD10->scale(xSec / hArea);
    } 

  }
  
}
