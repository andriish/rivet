// -*- C++ -*-
#include "Rivet/Tools/CentralityHistogram.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Analysis.hh"

namespace Rivet {

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
  double cest = centin >= 0.0? centin:
    applyProjection<CentralityEstimator>(ev, _centralityEstimatorName).estimate();
  _weightsum += ev.weight();

  // If estimator is negative, something has gone wrong.
  if ( cest < 0.0 ) {
    _currentHist = Histo1DPtr();
    return;
  }

  // If the current value of the centrality estimator is already in
  // use, just return the corresponding histogram.
  FlexiBinSet::iterator sit = _findBin(cest);
  if ( sit != _flexiBins.end() ) {
    sit->_weightsum += ev.weight();
    return;
  }

  if ( _flexiBins.size() >= _originalHists.size()*_overSamplingFactor )
    _currentHist = _purgeSamplers(); // Purge the histogram furthest
                                     // away from a centrality limit.
  else
    _currentHist = _newSamplerHist(); // Create a new sampler
                                      // histogram

  _flexiBins.insert(FlexiBin(_currentHist, cest, ev.weight()));

}

Histo1DPtr CentralityHistogram::_newSamplerHist() {
  // This is an error;
  if ( _originalHists.empty() ) return Histo1DPtr();

  Histo1DPtr hist(_originalHists.begin()->second->newclone());
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
  double totweight = _weightsum*100.0;
  double accweight = 0.0;
  
  for ( FlexiBinSet::iterator curr = _flexiBins.begin();
        curr != _flexiBins.end(); ++curr ) {
    accweight += 100.0*curr->_weightsum;
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
      // to 100, for all the rest we choose the middle bin.
      if ( *cit0 == 0.0 || *citn == 100.0 || counter%2 ) ++idist;
      if ( *citn == 100.0 || counter%2 ) ++midit;
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
  for ( FlexiBinSet::iterator fb = _flexiBins.begin();
        fb != _flexiBins.end(); ++fb ) {
    double chi = clo + 100.0*fb->_weightsum/_weightsum;
    for ( OriginalMap::iterator oit = _originalHists.begin();
          oit != _originalHists.end(); ++oit ) {
      double olo = oit->first.first;
      double ohi = oit->first.second;
      if ( clo > ohi || chi <= olo ) continue;
      // If we only have partial overlap we need to scale
      double lo = max(olo, clo);
      double hi = min(ohi, chi);
      Histo1D h = fb->_hist->clone();
      h.scaleW((hi - lo)/(chi - clo));
      *(oit->second) += h;
      _histWeightSum[oit->second] += fb->_weightsum*(hi - lo)/(chi - clo);
    }
  }
  _flexiBins.clear();
}

void CentralityHistogram::normalizePerEvent() {
  for ( OriginalMap::iterator oit = _originalHists.begin();
        oit != _originalHists.end(); ++oit ) {
    double sumw = eventWeightSum(oit->second);
    if ( sumw > 0.0 ) oit->second->normalize(1.0/sumw);
  }
}

}
