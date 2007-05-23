// -*- C++ -*-
#include "Rivet/Analysis/HepEx9506012.hh"
#include "Rivet/RivetCLHEP.hh"
#include "Rivet/RivetAIDA.hh"

using namespace Rivet;
using namespace std;
using namespace CLHEP;


/// @todo Declare via Cuts mechanism?
const double HepEx9506012::xmin = -6.0;
const double HepEx9506012::xmax = 6.0;


void HepEx9506012::init() {
  hEtFlow = vector<AIDA::IHistogram1D *>(nbin);
  hEtFlowStat = vector<AIDA::IHistogram1D *>(nbin);
  nev = vector<double>(nbin);
  for ( int i = 0; i < nbin; ++i ) {
    string I(1, char('1' + i));
    hEtFlow[i] = bookHistogram1D(I, "dEt/d[c] CMS bin=" + I, nb, xmin, xmax);
    hEtFlowStat[i] = bookHistogram1D(I, "stat dEt/d[c] CMS bin=1" + I, nb, xmin, xmax);
  }
  hAvEt = bookHistogram1D("21tmp", "<Et> vs kin. bin", nbin, 1.0, 10.0);
  hAvX = bookHistogram1D("22tmp", "<x>  vs kin. bin", nbin, 1.0, 10.0);
  hAvQ2 = bookHistogram1D("23tmp", "<Q2> vs kin. bin", nbin, 1.0, 10.0);
  hN = bookHistogram1D("24", "# events vs kin. bin", nbin, 1.0, 10.0);
}


int HepEx9506012::getbin(const DISKinematics & dk) {
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


void HepEx9506012::analyze(const Event & event) {
  const FinalStateHCM& fs = event.applyProjection(fsproj);
  const DISKinematics& dk = event.applyProjection(diskin);
  const CentralEtHCM y1 = event.applyProjection(y1hcm);
  
  int ibin = getbin(dk);
  if (ibin < 0) return;
  
  for ( int i = 0, N = fs.particles().size(); i < N; ++i ) {
    double rap = fs.particles()[i].getMomentum().rapidity();
    double et = fs.particles()[i].getMomentum().et();
    hEtFlow[ibin]->fill(rap, et*event.weight()/GeV);
    hEtFlowStat[ibin]->fill(rap, et*event.weight()/GeV);
  }
  
  nev[ibin] += event.weight();  
  hAvEt->fill(ibin + 1.5, y1.sumEt()*event.weight()/GeV);
  hAvX->fill(ibin + 1.5, dk.x()*event.weight());
  hAvQ2->fill(ibin + 1.5, dk.Q2()*event.weight()/(GeV*GeV));
  hN->fill(ibin + 1.5, event.weight());
}


void HepEx9506012::finalize() {
  for ( int ibin = 0; ibin < nbin; ++ibin ) {
    hEtFlow[ibin]->scale(1.0/(nev[ibin]*double(nb)/(xmax-xmin)));
    hEtFlowStat[ibin]->scale(1.0/(nev[ibin]*double(nb)/(xmax-xmin)));
  }
  AIDA::IHistogram1D* h = 0;
  h = histogramFactory().divide("/HepEx9506012/21", *hAvEt, *hN);
  h->setTitle(hAvEt->title());
  histogramFactory().destroy(hAvEt);
  h = histogramFactory().divide("/HepEx9506012/22", *hAvX, *hN);
  h->setTitle(hAvX->title());
  histogramFactory().destroy(hAvX);
  h = histogramFactory().divide("/HepEx9506012/23", *hAvQ2, *hN);
  h->setTitle(hAvQ2->title());
  histogramFactory().destroy(hAvQ2);
}
