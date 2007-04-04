// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/HepEx0409040.hh"
using namespace Rivet;

#include "AIDA/IHistogram1D.h"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////


// Book histograms
void HepEx0409040::init() {
  double bins[] = {64., 80., 88., 96., 100., 102., 104., 106., 108., 110.,
		    112., 114., 116., 118., 120., 122., 124., 126., 128.};
  histJetAzimuthpTmax75_100 = bookHistogram1D("JetAzimuthpTmax75_100", "Jet Jet azimuthal angle, pTmax=75..100", 18, 64., 128.);
  histJetAzimuthpTmax100_130 = bookHistogram1D("JetAzimuthpTmax100_130", "Jet Jet azimuthal angle, pTmax=100..130", 18, 64., 128.);
  histJetAzimuthpTmax130_180 = bookHistogram1D("JetAzimuthpTmax130_180", "Jet Jet azimuthal angle, pTmax=130..180", 18, 64., 128.);
  histJetAzimuthpTmax180_ = bookHistogram1D("JetAzimuthpTmax180_", "Jet Jet azimuthal angle, pTmax>180", 18, 64., 128.);
  /*
  histJetAzimuthpTmax75_100 = bookHistogram1D("JetAzimuthpTmax75_100", "Jet Jet azimuthal angle, pTmax=75..100", 18, bins);
  histJetAzimuthpTmax100_130 = bookHistogram1D("JetAzimuthpTmax100_130", "Jet Jet azimuthal angle, pTmax=100..130", 18, bins);
  histJetAzimuthpTmax130_180 = bookHistogram1D("JetAzimuthpTmax130_180", "Jet Jet azimuthal angle, pTmax=130..180", 18, bins);
  histJetAzimuthpTmax180_ = bookHistogram1D("JetAzimuthpTmax180_", "Jet Jet azimuthal angle, pTmax>180", 18, bins);
  */
}


// Do the analysis
void HepEx0409040::analyze(const Event & event) {
  Log& log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;
  
  // Analyse and print some info
  
  const D0RunIIconeJets& jetpro = event.applyProjection(conejets);
  
  int nj = jetpro.getNJets();
  
  //log << Log::INFO << "Jet multiplicity before any pT cut = " << nj << endl;
  
   
  //check fabs(z-vertex) < 50 cm 
  //const PVertex& PV = event.applyProjection(p_vertex); //segmentation violation
  //either the HepMC event record is not filled properly or the F77-Wrapper functions are faulty
  //if (fabs(PV().position().z())< 500.) { //assummed to be in mm

    
    // Fill histograms
    //for (std::vector<KtJet::KtLorentzVector>::iterator j = jetList.begin(); j != jetList.end(); j++) {
    //std::list<HepEntity>::iterator jetpTmax = jetpro.jets->begin(), jet2ndpTmax = jetpTmax;
    std::list<HepEntity>::iterator jetpTmax=0, jet2ndpTmax=0;
    //cout << "jetlist size = " << jetpro.jets->size() << endl;
    
    int Njet=0;
    for (std::list<HepEntity>::iterator jt = jetpro.jets->begin(); jt != jetpro.jets->end(); jt++) {
      //cout << "list item pT = " << jt->pT() << " E=" << jt->E << " pz=" << jt->pz << endl;
      if (jt->pT()>40.) {
	Njet++;
	//cout << "jet pT=" << jt->pT() << " y=" << jt->y() << " phi=" << jt->phi() << endl; 
      }
      if (jetpTmax==0 || jt->pT() > jetpTmax->pT()) {
	jet2ndpTmax = jetpTmax;
	jetpTmax = jt;
      }
      else if (jet2ndpTmax==0 || jt->pT() > jet2ndpTmax->pT()) {
	jet2ndpTmax = jt;
      }
    }
    
    //if (Njet>=2) log << Log::INFO << "Jet multiplicity after pT>40GeV cut = " << Njet << endl;
    
    /*
    if (jetpro.jets->size()>=2) {
      cout << "1st jet: E=" << jetpTmax->E << "  pz=" << jetpTmax->pz 
	   << "  2nd jet E=" << jet2ndpTmax->E << "  pz=" << jetpTmax->pz << endl;
      //cout << "1st jet: pT=" << jetpTmax->pT() << "  y=" << jetpTmax->y() 
      //   << "  2nd jet pT=" << jet2ndpTmax->pT() << "  y=" << jetpTmax->y() << endl;
    }
    */

    if (jetpro.jets->size()>=2 && jet2ndpTmax->pT() > 40.) {
      if (fabs(jetpTmax->y())<0.5 && fabs(jet2ndpTmax->y())<0.5) {
	//cout << "jet eta and pT requirements fulfilled" << endl;
	double etaMax = 3.0; //D0 calorimeter boundary
	bool addMuons = false; //Muons pass calorimeter almost without energy loss
	p_calmet.initialize(etaMax, addMuons);
	const CalMET& CaloMissEt = event.applyProjection(p_calmet);
	//cout << "CaloMissEt.MET()=" << CaloMissEt.MET() << endl;
	if (CaloMissEt.MET() < 0.7*jetpTmax->pT()) {
	  
	  double dphi = delta_phi(jetpTmax->phi(),jet2ndpTmax->phi());
	  dphi *= 128./PI; //publication histogramming choice
	  
	  if (jetpTmax->pT() > 75. && jetpTmax->pT() <= 100.)
	    histJetAzimuthpTmax75_100->fill(dphi, 1.0);
	  else if (jetpTmax->pT() > 100. && jetpTmax->pT() <= 130.)
	    histJetAzimuthpTmax100_130->fill(dphi, 1.0);
	  else if (jetpTmax->pT() > 130. && jetpTmax->pT() <= 180.)
	    histJetAzimuthpTmax130_180->fill(dphi, 1.0);
	  else if (jetpTmax->pT() > 180.)
	    histJetAzimuthpTmax180_->fill(dphi, 1.0);

	} //CalMET
      } //jets N, pT
    } //jets y (raqpidity) 
    
    //} //z-vertex
  

  // Finished...
  log << Log::DEBUG << "Finished analyzing" << endl;
}


// Finalize
void HepEx0409040::finalize() 
{ }


// Provide info object
RivetInfo HepEx0409040::getInfo() const {
  return Analysis::getInfo() + conejets.getInfo();
}
