// -*- C++ -*-
#ifndef RIVET_BINNEDHISTOGRAM_HH
#define RIVET_BINNEDHISTOGRAM_HH
#include "Rivet/Rivet.hh"

namespace Rivet {


  /**
   * BinnedHistogram contains a series of histograms of the same quantity
   * each in a different region of a second quantity.  For example, a
   * BinnedHistogram may contain histograms of the cross section differential
   * in PT in different eta regions.
   **/
  template<typename T> 
  class BinnedHistogram {
  public:
    
    /// Create a new empty BinnedHistogram
    BinnedHistogram() {
      return;
    }
 
    ///  Add a histogram in the region between @a binMin and @a binMax to this
    ///  set of BinnedHistograms.
    const BinnedHistogram<T>& addHistogram(const T& binMin,
                                           const T& binMax,
                                           AIDA::IHistogram1D* histo);
    
    /// Fill the histogram that lies in the same region as @a bin with the value
    /// @a val of weight @a weight.
    AIDA::IHistogram1D* fill(const T& bin,
                             const T& val,
                             double weight);
    
    const vector<AIDA::IHistogram1D*>& getHistograms() const { return _histos; }
    vector<AIDA::IHistogram1D*>& getHistograms() { return _histos; }
 

  private:
 
    map<T, AIDA::IHistogram1D*> _histosByUpperBound;
    map<T, AIDA::IHistogram1D*> _histosByLowerBound;
    vector<AIDA::IHistogram1D*> _histos;
 
  };


}

#endif
