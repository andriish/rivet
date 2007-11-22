// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Analyses/S3167097.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {

  const double S3167097::_xmin = -6.0;
  const double S3167097::_xmax = 6.0;


  void S3167097::init() {
    _hEtFlow = vector<AIDA::IHistogram1D *>(_nbin);
    _hEtFlowStat = vector<AIDA::IHistogram1D *>(_nbin);
    _nev = vector<double>(_nbin);
    /// @todo Automate this sort of thing so that the analysis code is more readable.
    for (int i = 0; i < _nbin; ++i) {
      string istr(1, char('1' + i));
      _hEtFlow[i] = bookHistogram1D(istr, "dEt/d[c] CMS bin=" + istr, _nb, _xmin, _xmax);
      _hEtFlowStat[i] = bookHistogram1D(istr, "stat dEt/d[c] CMS bin=1" + istr, _nb, _xmin, _xmax);
    }
    _hAvEt = bookHistogram1D("21tmp", "<Et> vs kin. bin", _nbin, 1.0, 10.0);
    _hAvX  = bookHistogram1D("22tmp", "<x>  vs kin. bin", _nbin, 1.0, 10.0);
    _hAvQ2 = bookHistogram1D("23tmp", "<Q2> vs kin. bin", _nbin, 1.0, 10.0);
    _hN    = bookHistogram1D("24", "# events vs kin. bin", _nbin, 1.0, 10.0);
  }
  

  int S3167097::getbin(const DISKinematics& dk) {
    const double GeV2 = GeV*GeV;

    if ( dk.Q2() > 5.0*GeV2 && dk.Q2() <= 10.0*GeV2 ) {
      if ( dk.x() > 0.0001 && dk.x() <= 0.0002 )
        return 0;
      else if ( dk.x() > 0.0002 && dk.x() <= 0.0005 && dk.Q2() > 6.0*GeV2 )
        return 1;
    }
    else if ( dk.Q2() > 10.0*GeV2 && dk.Q2() <= 20.0*GeV2 ){
      if ( dk.x() > 0.0002 && dk.x() <= 0.0005 )
        return 2;
      else if ( dk.x() > 0.0005 && dk.x() <= 0.0008 )
        return 3;
      else if ( dk.x() > 0.0008 && dk.x() <= 0.0015 )
        return 4;
      else if ( dk.x() > 0.0015 && dk.x() <= 0.0040 )
        return 5;
    }
    else if ( dk.Q2() > 20.0*GeV2 && dk.Q2() <= 50.0*GeV2 ){
      if ( dk.x() > 0.0005 && dk.x() <= 0.0014 )
        return 6;
      else if ( dk.x() > 0.0014 && dk.x() <= 0.0030 )
        return 7;
      else if ( dk.x() > 0.0030 && dk.x() <= 0.0100 )
        return 8;
    }
    return -1;
  }


  void S3167097::analyze(const Event& event) {
    const FinalStateHCM& fs = event.applyProjection(_fshcmproj);
    const DISKinematics& dk = event.applyProjection(_diskinproj);
    const CentralEtHCM y1 = event.applyProjection(_y1hcmproj);

    int ibin = getbin(dk);
    if (ibin < 0) return;

    const double weight = event.weight();

    for ( int i = 0, N = fs.particles().size(); i < N; ++i ) {
      const double rap = fs.particles()[i].getMomentum().rapidity();
      const double et = fs.particles()[i].getMomentum().Et();
      _hEtFlow[ibin]->fill(rap, weight * et/GeV);
      _hEtFlowStat[ibin]->fill(rap, weight * et/GeV);
    }

    _nev[ibin] += weight;
    _hAvEt->fill(ibin + 1.5, weight * y1.sumEt()/GeV);
    _hAvX->fill(ibin + 1.5, weight * dk.x());
    _hAvQ2->fill(ibin + 1.5, weight * dk.Q2()/(GeV*GeV));
    _hN->fill(ibin + 1.5, weight);
  }


  void S3167097::finalize() {
    for ( int ibin = 0; ibin < _nbin; ++ibin ) {
      _hEtFlow[ibin]->scale(1.0/(_nev[ibin]*double(_nb)/(_xmax-_xmin)));
      _hEtFlowStat[ibin]->scale(1.0/(_nev[ibin]*double(_nb)/(_xmax-_xmin)));
    }
    /// @todo Automate this sort of thing so that the analysis code is more readable.
    /// @todo Eliminate AIDA in favour of a histo interface that is "more C++", less factory-obsessive.
    AIDA::IHistogram1D* h = 0;
    h = histogramFactory().divide("/S3167097/21", *_hAvEt, *_hN);
    h->setTitle(_hAvEt->title());
    histogramFactory().destroy(_hAvEt);
    h = histogramFactory().divide("/S3167097/22", *_hAvX, *_hN);
    h->setTitle(_hAvX->title());
    histogramFactory().destroy(_hAvX);
    h = histogramFactory().divide("/S3167097/23", *_hAvQ2, *_hN);
    h->setTitle(_hAvQ2->title());
    histogramFactory().destroy(_hAvQ2);
  }

}
