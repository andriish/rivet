// -*- C++ -*-
#ifndef RIVET_CENTRALITYHISTOGRAM_HH
#define RIVET_CENTRALITYHISTOGRAM_HH
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Projections/CentralityEstimator.hh"
namespace Rivet {

class Analysis;


/**
 * CentralityHistogram contains a series of histograms of the same
 * quantity each in a different percentiles of a second quantity.
 * For example, a CentralityHistogram may contain histograms of the
 * cross section differential in \f$ p_T \f$ in different centrality
 * regions for heavy ion collisions.
 **/
class CentralityHistogram: public ProjectionApplier {
  public:

  /// Create a new empty CentralityHistogram. As an optional
  /// argument a CentralityProjector may be supplied. This will then
  /// be responsible for calculating a centrality estimate for each
  /// event.
  CentralityHistogram()
    : _currentHist(Histo1DPtr()), _currentCEst(-1.0),
      _overSamplingFactor(10), _weightsum(0.0) {
    _percentiles.insert(0.0);
    _percentiles.insert(1.0);
  }

  /// Set the centrality projection to be used. Note that this
  /// projection must have already been declared to Rivet.
  void setProjection(string pname) {
    _estimator = pname;
  }

  /// Return the class name.
  std::string name() const {
    return "Rivet::CentralityHistogram";
  }
  

  /// Add a histogram in the region between @a cmin and @a cmax to
  /// this set of CentralityHistograms. The range represent
  /// percentiles and must be between 0 and 100. No overlaping bins
  /// are allowed.

  /// Note that (cmin=0, cmax=5), means the five percent LEAST
  /// central events.
  void add(Histo1DPtr hist, double cmin, double cmax,
           double cestmin = -1.0, double cestmax = -1.0 );

  /// Setup the CentralityHistogram for the given event. Must be
  /// called for every event before any fill. Optionally an explicit
  /// value of the centrality estimator can be given, Otherwise a
  /// CentralityEstimator must have been provided beforehand. 
  void setup(const Event & ev, double cest = -1.0);

  /// Fill the histogram that lies in the same region as @a bin with the value
  /// @a x of weight @a weight.
  void fill(double x, double weight = 1.0) {
    if ( _currentHist ) _currentHist->fill(x, weight);
    for ( auto h : _ready )
      if ( h.second.inRange(_currentCEst) )
        h.second._hist->fill(x, weight);
  }

  /// Fill histo bin i with the given weight
  void fillBin(size_t i, double weight = 1.0) {
    if ( _currentHist ) _currentHist->fill(i, weight);
    for ( auto h : _ready )
      if ( h.second.inRange(_currentCEst) )
        h.second._hist->fill(i, weight);
  }

  /// At the end of the run, calculate the percentiles and fill the
  /// histograms provided with addHistogram(). This is typically
  /// called from the finalize method in a Analysis, but can also be
  /// called earlier in which case the the fill function can be
  /// continued to run as before.
  void finalize();

  /// Normalize each histogram to the sum of event weights in the
  /// corresponding centrality bin.
  void normalizePerEvent();

  /// Return a map bin edges of the centrality extimator indexed by
  /// the corresponing percentile.
  map<double,double> edges() const;

protected:

  /// Create a copy of the first original histogram to be used in
  /// the dynamical binning of centrality.
  Histo1DPtr _newSamplerHist();

  /// Take the sampler histogram which is furthest away from a
  /// centrality limit and merging it into the closest one and clear
  /// it and return it to be used for another dynamic bin.
  Histo1DPtr _purgeSamplers();
    

private:

  /// A flexible bin struct to be used to store temporary histograms.
  struct FlexiBin {

    /// Construct with an initial centrality estimate and an event
    /// weight.
    FlexiBin(Histo1DPtr hist, double cest = 0.0, double weight = 0.0)
      : _hist(hist), _cestLo(cest), _cestHi(cest), _weightsum(weight) {}

    /// Construct a temporary FlexiBin for finding a bin in a set.
    FlexiBin(double cest)
      : _cestLo(cest), _cestHi(cest), _weightsum(0.0) {}

    /// Merge in the contents of another FlexiBin into this. 
    void merge(const FlexiBin & fb) {
      _cestLo = min(_cestLo, fb._cestLo);
      _cestHi = max(_cestHi, fb._cestHi);
      _weightsum += fb._weightsum;
      *_hist += *fb._hist;
    }

    /// Comparisons for containers.
    bool operator< (const FlexiBin & fb) const {
      return ( _cestLo < fb._cestLo ||
               ( _cestLo == fb._cestLo && _cestLo < fb._cestHi ) );
    }

    /// The associated histogram.
    Histo1DPtr _hist;

    /// Current lower and upper edge of the centrality estimator for
    /// the fills in the associated histogram.
    double _cestLo, _cestHi;

    /// The sum of weights for all events entering the associated
    /// histogram.
    mutable double _weightsum;

  };

  struct HistBin {

    /// Construct a completely empty bin.
    HistBin()
      : _centLo(-1.0), _centHi(-1.0),
        _cestLo(-1.0), _cestHi(-1.0), _weightsum(0.0) {}

    /// Constructor taking a histogram and centrality interval as
    /// argument. Optionally the interval in the estimator can be
    /// given, in which case this histogram is considered to be
    /// "final".
    HistBin(Histo1DPtr h, double centLo, double centHi,
            double cestLo = -1.0, double cestHi = -1.0)
      : _hist(h), _centLo(centLo), _centHi(centHi),
        _cestLo(cestLo), _cestHi(cestHi), _weightsum(0.0) {}

    /// Return true if the given centrality estimate is in the range
    /// of this histogram.
    bool inRange(double cest) const {
      return _cestHi > 0 && _cestLo < cest && cest <= _cestHi;
    }

    /// Normalise the histogram to the tital cross section.
    void normalizePerEvent() {
      if ( _weightsum > 0.0 ) _hist->normalize(1.0/_weightsum);
    }

    /// The histogram.
    Histo1DPtr _hist;

    /// The range in centrality.
    double _centLo, _centHi;

    /// The corresponding range in the centrality estimator.
    double _cestLo, _cestHi;

    /// The sum of event weights for this bin;
    double _weightsum;

  };

  /// Convenient typedefs.
  typedef set<FlexiBin> FlexiBinSet;
  
  /// Find a bin corresponding to a given value of the centrality
  /// estimator.
  FlexiBinSet::iterator _findBin(double cest) const;

  /// The name of the CentralityEstimator projection to be used.
  string _estimator;

  /// The current temporary histogram selected for the centrality
  /// estimator calculated from the event presented in setup().
  Histo1DPtr _currentHist;

  /// THe current value of the centrality estimator.
  double _currentCEst;

  /// The oversampling of centrality bins. For each requested
  /// centrality bin this number of dynamic bins will be used.
  int _overSamplingFactor;

  /// The unfilled histograms where the esimator edges has not yet
  /// been determined.
  vector<HistBin> _unfilled;

  /// The dynamic bins for ranges of centrality estimators.
  FlexiBinSet _flexiBins;

  /// The sum of all event weights so far.
  double _weightsum;

  /// The sum of weights per original histogram.
  map<Histo1DPtr, double> _histWeightSum;

  /// Percentile limits.
  set<double> _percentiles;

  /// The filled histograms where the estimator edges has been determined.
  map<Histo1DPtr, HistBin> _ready;

};

}

#endif
