// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HZ95108 class.
//

#include "HZ95108.h"
#include "Rivet/Analysis/RivetHandler.h"
#include "Rivet/CLHEPWrap/SystemOfUnits.h"

using namespace Rivet;
using namespace std;

HZ95108::~HZ95108() {}

const double HZ95108::xmin = -6.0;
const double HZ95108::xmax = 6.0;

void HZ95108::init() {
  hEtFlow = vector<AIDA::IHistogram1D *>(nbin);
  hEtFlowStat = vector<AIDA::IHistogram1D *>(nbin);
  etmean = xmean = q2mean = nev = vector<double>(nbin);

  tree().mkdir("/HZ95108");
  for ( int i = 0; i < nbin; ++i ) {
    string I(1, char('1' + i));
    hEtFlow[i] =
      histogramFactory().createHistogram1D("/HZ95108/" + I, nb, xmin,xmax);
    hEtFlow[i]->setTitle("dEt/d[c] CMS bin=" +I);
    histogramFactory().createHistogram1D("/HZ95108/1" + I, nb, xmin,xmax);
    hEtFlow[i]->setTitle("stat dEt/d[c] CMS bin=1" + I);
  }
  hAvEt =
    histogramFactory().createHistogram1D("/HZ95108/21tmp", nbin, 1.0, 10.0);
  hAvX =
    histogramFactory().createHistogram1D("/HZ95108/22tmp", nbin, 1.0, 10.0);
  hAvX->setTitle("<x>  vs kin. bin");
  hAvQ2 =
    histogramFactory().createHistogram1D("/HZ95108/23tmp", nbin, 1.0, 10.0);
  hAvQ2->setTitle("<Q2> vs kin. bin");
  hN =
    histogramFactory().createHistogram1D("/HZ95108/24", nbin, 1.0, 10.0);
  hN->setTitle("# events vs kin. bin");
}

int HZ95108::getbin(const DISKinematics & dk) {
  const double GeV2 = CLHEP::GeV*CLHEP::GeV;

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

void HZ95108::analyze(const Event & event) {
  const FinalStateHCM & fs = event(fsproj);
  const DISKinematics & dk = event(diskin);

  int ibin = getbin(dk);
  if ( ibin < 0 ) return;

  double etcent = 0.0;

  for ( int i = 0, N = fs.particles().size(); i < N; ++i ) {
    double rap = fs.particles()[i].momentum.rapidity();
    double et = fs.particles()[i].momentum.et();
    hEtFlow[ibin]->fill(rap, et*event.weight());
    hEtFlowStat[ibin]->fill(rap, et*event.weight());

    if ( abs(rap) < 0.5 ) {
      etmean[ibin] += et;
      etcent += et;
    }
  }

  nev[ibin] += event.weight();
  xmean[ibin] += dk.x();
  q2mean[ibin] += dk.Q2();

  hAvEt->fill(ibin + 1.5, etcent*event.weight());
  hAvX->fill(ibin + 1.5, dk.x()*event.weight());
  hAvQ2->fill(ibin + 1.5, dk.Q2()*event.weight());
  hN->fill(ibin + 1.5, event.weight());

}

void HZ95108::finalize() {

  for ( int ibin = 0; ibin < nbin; ++ibin ) {
    hEtFlow[ibin]->scale(1.0/(nev[ibin]*double(nb)/(xmax-xmin)));
    hEtFlowStat[ibin]->scale(1.0/(nev[ibin]*double(nb)/(xmax-xmin)));
  }
  AIDA::IHistogram1D * h =
    histogramFactory().divide("/HZ95108/21", *hAvEt, *hN);
  h->setTitle(hAvEt->title());
  histogramFactory().destroy(hAvEt);
  h =
    histogramFactory().divide("/HZ95108/22", *hAvX, *hN);
  h->setTitle(hAvX->title());
  histogramFactory().destroy(hAvX);
  h =
    histogramFactory().divide("/HZ95108/23", *hAvQ2, *hN);
  h->setTitle(hAvQ2->title());
  histogramFactory().destroy(hAvQ2);

}

RivetInfo HZ95108::getInfo() const {
  return AnalysisBase::getInfo() + fsproj.getInfo();
}
