// -*- C++ -*-

#include "Rivet/Analyses/HepEx0701051.hh"

#include "Rivet/RivetAIDA.hh"

using namespace AIDA;

///////////////////////////////////////////////////////////////////////////////

namespace Rivet{
  
  const double HepEx0701051::_ktRParam = 0.7;
  const double HepEx0701051::_jetMinPT = 54.0;
  
  // Book histos and set counters for number of events passed in each one
  void HepEx0701051::init() {
    //@todo set _xSecTot from generator
    //just a guess at the generated cross section for now!
    _xSecTot = 600.0;
    _eventsTried = 0.0;
    _histos[0.1] = bookHistogram1D(1,1,1,"eta &lt; 0.1");
    _histos[0.7] = bookHistogram1D(2,1,1,"0.1 &lt; eta &lt; 0.7");
    _histos[1.1] = bookHistogram1D(3,1,1,"0.7 &lt; eta &lt; 1.1");
    _histos[1.6] = bookHistogram1D(4,1,1,"1.1 &lt; eta &lt; 1.6");
    _histos[2.1] = bookHistogram1D(5,1,1,"1.6 &lt; eta &lt; 2.1");
    
    for(map<double, IHistogram1D*>::iterator histIt = _histos.begin();
        histIt != _histos.end();
        ++histIt){
      _eventsPassed[histIt->second] = 0.0;
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////
  
  void HepEx0701051::analyze(const Event& event) {
    
    event.applyProjection(_ktproj);
    
    vector<KtJet::KtLorentzVector> jetList = _ktproj.getJets();
    
    set<IHistogram1D*> passed;
    
    const double weight = event.weight();
    
    for(vector<KtJet::KtLorentzVector>::iterator jet = jetList.begin();
        jet != jetList.end();
        ++jet){
      
      double pt = jet->perp();
      
      if(pt > _jetMinPT){
        map<double, IHistogram1D*>::iterator histIt = 
        _histos.upper_bound(jet->eta());
        
        if(histIt != _histos.end()){
          IHistogram1D* histo = histIt->second;
          histo->fill(pt, weight);
          if(histo->coordToIndex(pt) != IAxis::OVERFLOW_BIN){
            passed.insert(histo);
          }
        }
      }
    }
    
    //increment the event counters for each histogram
    _eventsTried += weight;
    
    for(set<IHistogram1D*>::iterator histIt = passed.begin();
        histIt != passed.end();
        ++histIt){
      
      _eventsPassed[*histIt] += weight;
    }
    return;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  
  void HepEx0701051::finalize() {
    
    //normalise histograms to cross section
    
    //double xSecPerEvent = _xSecTot / _eventsTried;
    double xSecPerEvent = crossSection() / _eventsTried;
    //HepData data is in nb, crossSection returns pb.
    xSecPerEvent = 0.001 * xSecPerEvent; 
    
    for(map<IHistogram1D*, double>::iterator histIt = _eventsPassed.begin();
        histIt != _eventsPassed.end();
        ++histIt){
      IHistogram1D* hist = histIt->first;
      double xSec = xSecPerEvent * histIt->second;
      int nBins = hist->axis().bins();
      double hArea = 0.0;
      for(int iBin = 0; iBin != nBins; ++iBin){
        hArea += hist->binHeight(iBin) * hist->axis().binWidth(iBin);
      }
      hist->scale(xSec / hArea);
    }
    return;
  }
  
}
///////////////////////////////////////////////////////////////////////////////
