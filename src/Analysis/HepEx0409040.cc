// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/HepEx0409040.hh"
using namespace Rivet;

#include "AIDA/IHistogram1D.h"
#include "AIDA/IAxis.h"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////


// Book histograms
void HepEx0409040::init() {

  vector<double> bins_pTmax75_100(19);
  bins_pTmax75_100[0] = 64.;
  bins_pTmax75_100[1] = 80.;
  bins_pTmax75_100[2] = 88.;
  bins_pTmax75_100[3] = 96.;
  bins_pTmax75_100[4] = 100.;
  bins_pTmax75_100[5] = 102.;
  bins_pTmax75_100[6] = 104.;
  bins_pTmax75_100[7] = 106.;
  bins_pTmax75_100[8] = 108.;
  bins_pTmax75_100[9] = 110.;
  bins_pTmax75_100[10] = 112.;
  bins_pTmax75_100[11] = 114.;
  bins_pTmax75_100[12] = 116.;
  bins_pTmax75_100[13] = 118.;
  bins_pTmax75_100[14] = 120.;
  bins_pTmax75_100[15] = 122.;
  bins_pTmax75_100[16] = 124.;
  bins_pTmax75_100[17] = 126.;
  bins_pTmax75_100[18] = 128.;

  vector<double> bins_pTmax100_130(23);
  bins_pTmax100_130[0] = 64.;
  bins_pTmax100_130[1] = 72.;
  bins_pTmax100_130[2] = 80.;
  bins_pTmax100_130[3] = 84.;
  bins_pTmax100_130[4] = 88.;
  bins_pTmax100_130[5] = 92.;
  bins_pTmax100_130[6] = 94.;
  bins_pTmax100_130[7] = 96.;
  bins_pTmax100_130[8] = 98.;
  bins_pTmax100_130[9] = 100.;
  bins_pTmax100_130[10] = 102.;
  bins_pTmax100_130[11] = 104.;
  bins_pTmax100_130[12] = 106.;
  bins_pTmax100_130[13] = 108.;
  bins_pTmax100_130[14] = 110.;
  bins_pTmax100_130[15] = 112.;
  bins_pTmax100_130[16] = 114.;
  bins_pTmax100_130[17] = 116.;
  bins_pTmax100_130[18] = 118.;
  bins_pTmax100_130[19] = 120.;
  bins_pTmax100_130[20] = 122.;
  bins_pTmax100_130[21] = 124.;
  bins_pTmax100_130[22] = 126.;
  bins_pTmax100_130[23] = 128.;

  vector<double> bins_pTmax130_180(27);
  bins_pTmax130_180[0] = 64.;
  bins_pTmax130_180[1] = 80.;
  bins_pTmax130_180[2] = 88.;
  bins_pTmax130_180[3] = 92.;
  bins_pTmax130_180[4] = 96.;
  bins_pTmax130_180[5] = 98.;
  bins_pTmax130_180[6] = 100.;
  bins_pTmax130_180[7] = 102.;
  bins_pTmax130_180[8] = 104.;
  bins_pTmax130_180[9] = 106.;
  bins_pTmax130_180[10] = 108.;
  bins_pTmax130_180[11] = 110.;
  bins_pTmax130_180[12] = 112.;
  bins_pTmax130_180[13] = 114.;
  bins_pTmax130_180[14] = 115.;
  bins_pTmax130_180[15] = 116.;
  bins_pTmax130_180[16] = 117.;
  bins_pTmax130_180[17] = 118.;
  bins_pTmax130_180[18] = 119.;
  bins_pTmax130_180[19] = 120.;
  bins_pTmax130_180[20] = 121.;
  bins_pTmax130_180[21] = 122.;
  bins_pTmax130_180[22] = 123.;
  bins_pTmax130_180[23] = 124.;
  bins_pTmax130_180[24] = 125.;
  bins_pTmax130_180[25] = 126.;
  bins_pTmax130_180[26] = 127.;
  bins_pTmax130_180[27] = 128.;

  vector<double> bins_pTmax180_(26);
  bins_pTmax180_[0] = 64.;
  bins_pTmax180_[1] = 80.;
  bins_pTmax180_[2] = 88.;
  bins_pTmax180_[3] = 92.;
  bins_pTmax180_[4] = 96.;
  bins_pTmax180_[5] = 98.;
  bins_pTmax180_[6] = 100.;
  bins_pTmax180_[7] = 102.;
  bins_pTmax180_[8] = 104.;
  bins_pTmax180_[9] = 106.;
  bins_pTmax180_[10] = 108.;
  bins_pTmax180_[11] = 110.;
  bins_pTmax180_[12] = 112.;
  bins_pTmax180_[13] = 114.;
  bins_pTmax180_[14] = 115.;
  bins_pTmax180_[15] = 116.;
  bins_pTmax180_[16] = 117.;
  bins_pTmax180_[17] = 118.;
  bins_pTmax180_[18] = 119.;
  bins_pTmax180_[19] = 120.;
  bins_pTmax180_[20] = 121.;
  bins_pTmax180_[21] = 122.;
  bins_pTmax180_[22] = 123.;
  bins_pTmax180_[23] = 124.;
  bins_pTmax180_[24] = 125.;
  bins_pTmax180_[25] = 126.;
  bins_pTmax180_[26] = 127.;



  /*  
  histJetAzimuthpTmax75_100 = bookHistogram1D("JetAzimuthpTmax75_100", "Jet Jet azimuthal angle, pTmax=75..100", 18, 64., 128.);
  histJetAzimuthpTmax100_130 = bookHistogram1D("JetAzimuthpTmax100_130", "Jet Jet azimuthal angle, pTmax=100..130", 18, 64., 128.);
  histJetAzimuthpTmax130_180 = bookHistogram1D("JetAzimuthpTmax130_180", "Jet Jet azimuthal angle, pTmax=130..180", 18, 64., 128.);
  histJetAzimuthpTmax180_ = bookHistogram1D("JetAzimuthpTmax180_", "Jet Jet azimuthal angle, pTmax>180", 18, 64., 128.);
  */
  histJetAzimuthpTmax75_100 = bookHistogram1D("JetAzimuthpTmax75_100", "Jet Jet azimuthal angle, pTmax=75..100", bins_pTmax75_100);
  histJetAzimuthpTmax100_130 = bookHistogram1D("JetAzimuthpTmax100_130", "Jet Jet azimuthal angle, pTmax=100..130", bins_pTmax100_130);
  histJetAzimuthpTmax130_180 = bookHistogram1D("JetAzimuthpTmax130_180", "Jet Jet azimuthal angle, pTmax=130..180", bins_pTmax130_180);
  histJetAzimuthpTmax180_ = bookHistogram1D("JetAzimuthpTmax180_", "Jet Jet azimuthal angle, pTmax>180", bins_pTmax180_);

}


// Do the analysis
void HepEx0409040::analyze(const Event & event) {
  Log& log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;
  //cout << "Start analyzing" << endl; //ls

  // Analyse and print some info
  
  const D0RunIIconeJets& jetpro = event.applyProjection(conejets);
  
  //int nj = jetpro.getNJets(); //ls
  //log << Log::INFO << "Jet multiplicity before any pT cut = " << nj << endl; //ls
  //cout << "Jet multiplicity before any pT cut = " << nj << endl; //ls
   
  //check fabs(z-vertex) < 50 cm 
  const PVertex& PV = event.applyProjection(p_vertex); //segmentation violation
  //either the HepMC event record is not filled properly or the F77-Wrapper functions are faulty
  if (fabs(PV().position().z())< 500.) { //assummed to be in mm, PYTHIA convention

    
    // Fill histograms
    std::list<HepEntity>::iterator jetpTmax = jetpro.jets->end(),
      jet2ndpTmax = jetpro.jets->end();
    //cout << "jetlist size = " << jetpro.jets->size() << endl;
    
    int Njet=0;
    for (std::list<HepEntity>::iterator jt = jetpro.jets->begin();
	 jt != jetpro.jets->end(); jt++) {
      //cout << "list item pT = " << jt->pT() << " E=" << jt->E << " pz=" << jt->pz << endl;
      if (jt->pT()>40.) {
	Njet++;
	//cout << "jet pT=" << jt->pT() << " y=" << jt->y() << " phi=" << jt->phi() << endl; 
      }
      if (jetpTmax == jetpro.jets->end() || jt->pT() > jetpTmax->pT()) {
	jet2ndpTmax = jetpTmax;
	jetpTmax = jt;
      }
      else if (jet2ndpTmax == jetpro.jets->end() ||
	       jt->pT() > jet2ndpTmax->pT()) {
	jet2ndpTmax = jt;
      }
    }
       
    //if (Njet>=2) {
    //log << Log::INFO << "Jet multiplicity after pT>40GeV cut = " << Njet << endl; //ls
    //cout << "Jet multiplicity after pT>40GeV cut = " << Njet << endl; //ls
    //}

    /*
    if (jetpro.jets->size()>=2) {
      cout << "1st jet: E=" << jetpTmax->E << "  pz=" << jetpTmax->pz << " pt=" << jetpTmax->pT() << endl;
      cout << "2nd jet E=" << jet2ndpTmax->E << "  pz=" << jet2ndpTmax->pz << " pt=" << jet2ndpTmax->pT() << endl;
      //cout << "1st jet: pT=" << jetpTmax->pT() << "  y=" << jetpTmax->y() 
      //   << "  2nd jet pT=" << jet2ndpTmax->pT() << "  y=" << jetpTmax->y() << endl;
    }
    */

    if (jetpro.jets->size()>=2 && jet2ndpTmax->pT() > 40.) {
      if (fabs(jetpTmax->y())<0.5 && fabs(jet2ndpTmax->y())<0.5) {
	//cout << "jet eta and pT requirements fulfilled" << endl; //ls
	double etaMax = 3.0; //D0 calorimeter boundary
	bool addMuons = false; //Muons pass calorimeter almost without energy loss
	p_calmet.initialize(etaMax, addMuons);
	const CalMET& CaloMissEt = event.applyProjection(p_calmet);
	//cout << "CaloMissEt.MET()=" << CaloMissEt.MET() << endl; //ls
	if (CaloMissEt.MET() < 0.7*jetpTmax->pT()) {
	  
	  double dphi = delta_phi(jetpTmax->phi(),jet2ndpTmax->phi());
	  dphi *= 128./PI; //publication histogramming choice
	  
	  //cout << "Filling histograms now: dphi=" << dphi << endl;

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
    
  } //z-vertex
  

  // Finished...
  log << Log::DEBUG << "Finished analyzing" << endl;
}


// Finalize
void HepEx0409040::finalize() { 

  Log& log = getLog();

  double area75_100 = 0;
  for (int i=0; i < histJetAzimuthpTmax75_100->axis().bins(); ++i) {
    area75_100 += histJetAzimuthpTmax75_100->binHeight(i) * 
      histJetAzimuthpTmax75_100->axis().binWidth(i);
  }
  log << Log::INFO << "Area under histJetAzimuthpTmax75_100 histogram: " << area75_100 << endl;

  double area100_130 = 0;
  for (int i=0; i < histJetAzimuthpTmax100_130->axis().bins(); ++i) {
    area100_130 += histJetAzimuthpTmax100_130->binHeight(i) * 
      histJetAzimuthpTmax100_130->axis().binWidth(i);
  }
  log << Log::INFO << "Area under histJetAzimuthpTmax100_130 histogram: " << area100_130 << endl;

  double area130_180 = 0;
  for (int i=0; i < histJetAzimuthpTmax130_180->axis().bins(); ++i) {
    area130_180 += histJetAzimuthpTmax130_180->binHeight(i) * 
      histJetAzimuthpTmax130_180->axis().binWidth(i);
  }
  log << Log::INFO << "Area under histJetAzimuthpTmax130_180 histogram: " << area130_180 << endl;

  double area180_ = 0;
  for (int i=0; i < histJetAzimuthpTmax180_->axis().bins(); ++i) {
    area180_ += histJetAzimuthpTmax180_->binHeight(i) * 
      histJetAzimuthpTmax180_->axis().binWidth(i);
  }
  log << Log::INFO << "Area under histJetAzimuthpTmax180_ histogram: " << area180_ << endl;



  // //Normalize the histogram areas to 1
  //Normalize to cross section (= sum of all data bins per pT histogram)
  if (histJetAzimuthpTmax75_100->sumBinHeights()!=0){
    //histJetAzimuthpTmax75_100->scale(1/histJetAzimuthpTmax75_100->sumBinHeights() );
    histJetAzimuthpTmax75_100->scale(19.9778/histJetAzimuthpTmax75_100->sumBinHeights() );
  }
  if (histJetAzimuthpTmax100_130->sumBinHeights()!=0){
    //histJetAzimuthpTmax100_130->scale(1/histJetAzimuthpTmax100_130->sumBinHeights() );
    histJetAzimuthpTmax100_130->scale(20.2871/histJetAzimuthpTmax100_130->sumBinHeights() );
  }
  if (histJetAzimuthpTmax130_180->sumBinHeights()!=0){
    //histJetAzimuthpTmax130_180->scale(1/histJetAzimuthpTmax130_180->sumBinHeights() );
    histJetAzimuthpTmax130_180->scale(38.1651/histJetAzimuthpTmax130_180->sumBinHeights() );
  }
  if (histJetAzimuthpTmax180_->sumBinHeights()!=0){
    //histJetAzimuthpTmax180_->scale(1/histJetAzimuthpTmax180_->sumBinHeights() );
    histJetAzimuthpTmax180_->scale(38.86794/histJetAzimuthpTmax180_->sumBinHeights() );
  }


  /*
  //histSphericity_->scale(1/histSphericity_->sumBinHeights() );
  histSphericity_->scale(1/area);
  histPlanarity_->scale(1/histPlanarity_->sumBinHeights() );
  histAplanarity_->scale(1/histAplanarity_->sumBinHeights() );
  */

  log << Log::INFO << "Sum of histJetAzimuthpTmax75_100 bin heights after normalization: "
      << histJetAzimuthpTmax75_100->sumBinHeights() << endl;
  log << Log::INFO << "Sum of histJetAzimuthpTmax100_130 bin heights after normalization: "
      << histJetAzimuthpTmax100_130->sumBinHeights() << endl;
  log << Log::INFO << "Sum of histJetAzimuthpTmax130_180 bin heights after normalization: "
      << histJetAzimuthpTmax130_180->sumBinHeights() << endl;
  log << Log::INFO << "Sum of histJetAzimuthpTmax180_ bin heights after normalization: "
      << histJetAzimuthpTmax180_->sumBinHeights() << endl;


}


// // Provide info object
// RivetInfo HepEx0409040::getInfo() const {
//   return Analysis::getInfo() + conejets.getInfo();
// }
