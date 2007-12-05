// -*- C++ -*-

#include "Rivet/Analyses/CDF_2007_S7057202.hh"
#include "Rivet/RivetAIDA.hh"
using namespace AIDA;


namespace Rivet {
  
  // Book histos and set counters for number of events passed in each one
  void CDF_2007_S7057202::init() {

    _eventsTried = 0.0;
    /// @todo Indexing a map by double is a bad idea...
    /// @todo XML escaping should be done by LWH
    _histos[0.1] = bookHistogram1D(1,1,1,"eta &lt; 0.1");
    _histos[0.7] = bookHistogram1D(2,1,1,"0.1 &lt; eta &lt; 0.7");
    _histos[1.1] = bookHistogram1D(3,1,1,"0.7 &lt; eta &lt; 1.1");
    _histos[1.6] = bookHistogram1D(4,1,1,"1.1 &lt; eta &lt; 1.6");
    _histos[2.1] = bookHistogram1D(5,1,1,"1.6 &lt; eta &lt; 2.1");
    
    for (map<double, IHistogram1D*>::iterator histIt = _histos.begin();
        histIt != _histos.end(); ++histIt) {
      _eventsPassed[histIt->second] = 0.0;
    }
  }
  

  
  void CDF_2007_S7057202::analyze(const Event& event) {
    const double weight = event.weight();
    event.applyProjection(_jetproj);
    typedef vector<fastjet::PseudoJet> Jets;
    Jets jetList = _jetproj.getJets();

    /// @todo This is an immensely fiddly way to fill the appropriate histo...

    set<IHistogram1D*> passed;
    for (Jets::const_iterator jet = jetList.begin(); jet != jetList.end(); ++jet) {
      const double pt = jet->perp();
      
      // If pass the pT cut and in the eta range of a histo...
      if (pt > _jetMinPT) {
        map<double, IHistogram1D*>::iterator histIt = _histos.upper_bound(jet->eta());
        if (histIt != _histos.end()) {
          IHistogram1D* histo = histIt->second;
          // ...fill the histo
          histo->fill(pt, weight);
          if (histo->coordToIndex(pt) != IAxis::OVERFLOW_BIN) {
            passed.insert(histo);
          }
        }
      }
    }
    
    // Increment the event counters for each histogram
    _eventsTried += weight;
    for (set<IHistogram1D*>::iterator histIt = passed.begin(); histIt != passed.end(); ++histIt){
      _eventsPassed[*histIt] += weight;
    }
  }
  


  // Normalise histograms to cross-section  
  void CDF_2007_S7057202::finalize() {
    double xSecPerEvent = crossSection() / _eventsTried;
    /// HepData data is in nb, crossSection returns pb.
    xSecPerEvent = 0.001 * xSecPerEvent; 
    
    for (map<IHistogram1D*, double>::iterator histIt = _eventsPassed.begin();
        histIt != _eventsPassed.end(); ++histIt) {
      IHistogram1D* hist = histIt->first;
      double xSec = xSecPerEvent * histIt->second;
      /// @todo Can't the normalize() function be used here?
      size_t nBins = hist->axis().bins();
      double hArea = 0.0;
      for (size_t iBin = 0; iBin != nBins; ++iBin) {
        hArea += hist->binHeight(iBin) * hist->axis().binWidth(iBin);
      }
      hist->scale(xSec / hArea);
    }
  }

  
}
