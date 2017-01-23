// -*- C++ -*-
#include "Rivet/Tools/CentralityHistogram.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Analysis.hh"

namespace Rivet {

void CentralityHistogram::add(Histo1DPtr hist, double cmin, double cmax,
                              double cestmin, double cestmax) {
  if ( cestmax < 0.0 )
    _unfilled.push_back(HistBin(hist, cmin/100.0, cmax/100.0));
  else 
    _ready[hist] = HistBin(hist, cmin/100.0, cmax/100.0, cestmin, cestmax);
  _percentiles.insert(cmin/100.0);
  _percentiles.insert(cmax/100.0);
}


CentralityHistogram::FlexiBinSet::iterator
CentralityHistogram::_findBin(double cest) const {
  FlexiBinSet::iterator it = _flexiBins.upper_bound(FlexiBin(cest));
  if ( it != _flexiBins.end() && it != _flexiBins.begin() ) {
    --it;
    if ( it->_cestLo == cest ||
         ( it->_cestLo < cest && cest < it->_cestHi ) )
      return it;
  }
  return _flexiBins.end();
}

void CentralityHistogram::setup(const Event & ev, double centin) {

  // Get the estimator value from the associated projection.
  _currentCEst = centin >= 0.0? centin:
    applyProjection<CentralityEstimator>(ev, _estimator).estimate();
  _weightsum += ev.weight();
  for ( auto h : _ready )
    if ( h.second.inRange(_currentCEst) )
      h.second._weightsum += ev.weight();

  // If estimator is negative, something has gone wrong.
  if ( _currentCEst < 0.0 ) {
    _currentHist = Histo1DPtr();
    return;
  }

  // If the current value of the centrality estimator is already in
  // use, just return the corresponding histogram.
  FlexiBinSet::iterator sit = _findBin(_currentCEst);
  if ( sit != _flexiBins.end() ) {
    sit->_weightsum += ev.weight();
    return;
  }

  if ( _flexiBins.size() >= _unfilled.size()*_overSamplingFactor )
    _currentHist = _purgeSamplers(); // Purge the histogram furthest
                                     // away from a centrality limit.
  else
    _currentHist = _newSamplerHist(); // Create a new sampler
                                      // histogram

  _flexiBins.insert(FlexiBin(_currentHist, _currentCEst, ev.weight()));

}

Histo1DPtr CentralityHistogram::_newSamplerHist() {
  // This is an error;
  if ( _unfilled.empty() ) return Histo1DPtr();

  Histo1DPtr hist(_unfilled.begin()->_hist->newclone());
  hist->reset();
  return hist;

}

Histo1DPtr CentralityHistogram::_purgeSamplers() {

  set<double>::iterator citn = _percentiles.begin();
  set<double>::iterator cit0 = citn++;

  FlexiBinSet::iterator selectit = _flexiBins.begin();
  int selectdist = 0;
  int counter = 0;

  FlexiBinSet::iterator midit = selectit;
  int idist = 0;
  double totweight = _weightsum;
  double accweight = 0.0;
  
  for ( auto curr = _flexiBins.begin(); curr != _flexiBins.end(); ++curr ) {
    accweight += curr->_weightsum;
    if ( accweight > (*citn)*totweight ) {
      if ( idist > selectdist ) {
        selectdist = idist;
        selectit = midit;
      }
      idist = 0;
      midit = curr;
      cit0 = citn++;
      counter = 0;
    } else {
       ++counter;
      // Note that for the first centrality bin we always merge the
      // bin closest to 0, and for the last bin always the one closest
      // to 100%, for all the rest we choose the middle bin.
      if ( *cit0 == 0.0 || *citn == 1.0 || counter%2 ) ++idist;
      if ( *citn == 1.0 || counter%2 ) ++midit;
    }
  }

  FlexiBinSet::iterator mergeit = selectit;
  if ( selectit == _flexiBins.begin() ) ++selectit;
  else --mergeit;

  FlexiBin merged = *mergeit;
  FlexiBin selected = *selectit;
  merged.merge(*selectit);
  _flexiBins.erase(mergeit);
  _flexiBins.erase(selectit);
  _flexiBins.insert(merged);

  return selected._hist;
  
}

void CentralityHistogram::finalize() {

  // Take the contents of the dynamical binning and fill the original
  // histograms.

  double clo = 0.0;
  for ( const FlexiBin & fb : _flexiBins ) {
    double chi = clo + fb._weightsum/_weightsum;
    for ( HistBin & hbin : _unfilled ) {
      double olo = hbin._centLo;
      double ohi = hbin._centHi;
      if ( clo > ohi || chi <= olo ) continue;
      // If we only have partial overlap we need to scale
      double lo = max(olo, clo);
      double hi = min(ohi, chi);
      Histo1D h = fb._hist->clone();
      h.scaleW((hi - lo)/(chi - clo));
      *(hbin._hist) += h;
      hbin._weightsum += fb._weightsum*(hi - lo)/(chi - clo);
      if ( clo <= olo )
        hbin._cestLo =
          fb._cestLo + (fb._cestHi - fb._cestHi)*(olo - clo)/(chi - clo);
      if ( chi > ohi )
        hbin._cestHi =
          fb._cestLo + (fb._cestHi - fb._cestHi)*(ohi - clo)/(chi - clo);
    }
  }
  _flexiBins.clear();
  for ( HistBin & hbin : _unfilled ) _ready[hbin._hist] = hbin;
  _unfilled.clear();

}

void CentralityHistogram::normalizePerEvent() {
  for ( auto h : _ready ) h.second.normalizePerEvent();
}

map<double,double> CentralityHistogram::edges() const {
  map<double,double> ret;
  for ( auto hbin : _ready ) {
    ret[hbin.second._centLo] = hbin.second._cestLo;
    ret[hbin.second._centHi] = hbin.second._cestHi;
  }
  return ret;
}


}
