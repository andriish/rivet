// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/HepEx0409040.hh"
using namespace Rivet;

#include "Rivet/RivetAIDA.hh"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////


// Book histograms
void HepEx0409040::init() {

  //set in constructor (default=x-axis range in rad)
  //xscale = 1.;   // histo x-axis range in [rad*128/Pi]
  //xscale = PI/128.;   // histo x-axis range in [rad]


//   double bins_pTmax75_100_arr[] = {64., 80., 88., 96., 100., 102., 104., 106., 108., 110., 
// 				   112., 114., 116., 118., 120., 122., 124., 126., 128.}; 
//   vector<double> bins_pTmax75_100(19);
//   for (int i=0; i<=18; ++i) 
//     bins_pTmax75_100[i] = bins_pTmax75_100_arr[i]*xscale;

    
//   double bins_pTmax100_130_arr[] = {64., 72., 80., 84., 88., 92., 94., 96., 98., 
// 				    100., 102., 104., 106., 108., 110., 112., 114., 
// 				    116., 118., 120., 122., 124., 126., 128.}; 
//   vector<double> bins_pTmax100_130(24);
//   for (int i=0; i<=23; ++i)
//     bins_pTmax100_130[i] = bins_pTmax100_130_arr[i]*xscale;
  
  
//   double bins_pTmax130_180_arr[] = {64., 80., 88., 92., 96., 98., 100., 102., 104., 106., 
// 				    108., 110., 112., 114., 115., 116., 117., 118., 119.,
// 				    120., 121., 122., 123., 124., 125., 126., 127., 128.}; 
//   vector<double> bins_pTmax130_180(28);
//   for (int i=0; i<=27; ++i)
//     bins_pTmax130_180[i] = bins_pTmax130_180_arr[i]*xscale;
  
  
//   double bins_pTmax180_arr[] = {64., 80., 88., 92., 96., 98., 100., 102., 104., 106., 
// 				108., 110., 112., 114., 115., 116., 117., 118., 119., 
// 				120., 121., 122., 123., 124., 125., 126., 127.}; 
//   vector<double> bins_pTmax180_(27);
//   for (int i=0; i<=26; ++i) 
//     bins_pTmax180_[i] = bins_pTmax180_arr[i]*xscale;

  

//   histJetAzimuth_pTmax75_100 = bookHistogram1D("JetAzimuth_pTmax75_100", "Jet Jet azimuthal angle, pTmax=75..100", bins_pTmax75_100);
//   histJetAzimuth_pTmax100_130 = bookHistogram1D("JetAzimuth_pTmax100_130", "Jet Jet azimuthal angle, pTmax=100..130", bins_pTmax100_130);
//   histJetAzimuth_pTmax130_180 = bookHistogram1D("JetAzimuth_pTmax130_180", "Jet Jet azimuthal angle, pTmax=130..180", bins_pTmax130_180);
//   histJetAzimuth_pTmax180_ = bookHistogram1D("JetAzimuth_pTmax180_", "Jet Jet azimuthal angle, pTmax>180", bins_pTmax180_);


  // Use histogram auto-booking
  histJetAzimuth_pTmax75_100  = bookHistogram1D(1, 1, 1, "Jet Jet azimuthal angle, pTmax=75..100");
  histJetAzimuth_pTmax100_130 = bookHistogram1D(2, 1, 1, "Jet Jet azimuthal angle, pTmax=100..130");
  histJetAzimuth_pTmax130_180 = bookHistogram1D(3, 1, 1, "Jet Jet azimuthal angle, pTmax=130..180");
  histJetAzimuth_pTmax180_    = bookHistogram1D(4, 1, 1, "Jet Jet azimuthal angle, pTmax>180");
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
    list<HepEntity>::iterator jetpTmax = jetpro.jets->end();
    list<HepEntity>::iterator jet2ndpTmax = jetpro.jets->end();
    //cout << "jetlist size = " << jetpro.jets->size() << endl;
    
    int Njet=0;
    for (list<HepEntity>::iterator jt = jetpro.jets->begin(); jt != jetpro.jets->end(); ++jt) {
      //cout << "list item pT = " << jt->pT() << " E=" << jt->E << " pz=" << jt->pz << endl;
      if (jt->pT()>40.) ++Njet;
      //cout << "jet pT=" << jt->pT() << " y=" << jt->y() << " phi=" << jt->phi() << endl; 

      if (jetpTmax == jetpro.jets->end() || jt->pT() > jetpTmax->pT()) {
        jet2ndpTmax = jetpTmax;
        jetpTmax = jt;
      } else {
        if (jet2ndpTmax == jetpro.jets->end() || jt->pT() > jet2ndpTmax->pT()) {
          jet2ndpTmax = jt;
        }
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
          if (fabs(xscale-1.)<1.e-3) dphi /= xscale; //x-axis range [64,128]
          
          //cout << "Filling histograms now: dphi=" << dphi << endl;
          
          if (jetpTmax->pT() > 75. && jetpTmax->pT() <= 100.) 
            histJetAzimuth_pTmax75_100->fill(dphi, event.weight() );
          else if (jetpTmax->pT() > 100. && jetpTmax->pT() <= 130.)
            histJetAzimuth_pTmax100_130->fill(dphi, event.weight() );
          else if (jetpTmax->pT() > 130. && jetpTmax->pT() <= 180.)
            histJetAzimuth_pTmax130_180->fill(dphi, event.weight() );
          else if (jetpTmax->pT() > 180.)
            histJetAzimuth_pTmax180_->fill(dphi, event.weight() );
          
        } //CalMET
      } //jets N, pT
    } //jets y (raqpidity) 
    
  } //z-vertex
  
  
  // Finished
  log << Log::DEBUG << "Finished analyzing" << endl;
}


// Finalize
void HepEx0409040::finalize() { 
  Log& log = getLog();

//   double area75_100 = 0;
//   for (int i=0; i < histJetAzimuth_pTmax75_100->axis().bins(); ++i) {
//     area75_100 += histJetAzimuth_pTmax75_100->binHeight(i) * 
//       histJetAzimuth_pTmax75_100->axis().binWidth(i);
//   }
//   log << Log::INFO << "Area under histJetAzimuth_pTmax75_100 histogram: " << area75_100 << endl;

//   double area100_130 = 0;
//   for (int i=0; i < histJetAzimuth_pTmax100_130->axis().bins(); ++i) {
//     area100_130 += histJetAzimuth_pTmax100_130->binHeight(i) * 
//       histJetAzimuth_pTmax100_130->axis().binWidth(i);
//   }
//   log << Log::INFO << "Area under histJetAzimuth_pTmax100_130 histogram: " << area100_130 << endl;

//   double area130_180 = 0;
//   for (int i=0; i < histJetAzimuth_pTmax130_180->axis().bins(); ++i) {
//     area130_180 += histJetAzimuth_pTmax130_180->binHeight(i) * 
//       histJetAzimuth_pTmax130_180->axis().binWidth(i);
//   }
//   log << Log::INFO << "Area under histJetAzimuth_pTmax130_180 histogram: " << area130_180 << endl;

//   double area180_ = 0;
//   for (int i=0; i < histJetAzimuth_pTmax180_->axis().bins(); ++i) {
//     area180_ += histJetAzimuth_pTmax180_->binHeight(i) * 
//       histJetAzimuth_pTmax180_->axis().binWidth(i);
//   }
//   log << Log::INFO << "Area under histJetAzimuth_pTmax180_ histogram: " << area180_ << endl;



//   // //Normalize the histogram areas to 1
//   if (histJetAzimuth_pTmax75_100->sumBinHeights()!=0){
//     histJetAzimuth_pTmax75_100->scale(1./area75_100);
//   }
//   if (histJetAzimuth_pTmax100_130->sumBinHeights()!=0){
//     histJetAzimuth_pTmax100_130->scale(1./area100_130);
//   }
//   if (histJetAzimuth_pTmax130_180->sumBinHeights()!=0){
//     histJetAzimuth_pTmax130_180->scale(1./area130_180);
//   }
//   if (histJetAzimuth_pTmax180_->sumBinHeights()!=0){
//     histJetAzimuth_pTmax180_->scale(1./area180_); 
//   }

  // Normalize histograms to unit area
  normalize(histJetAzimuth_pTmax75_100);
  normalize(histJetAzimuth_pTmax100_130);
  normalize(histJetAzimuth_pTmax130_180);
  normalize(histJetAzimuth_pTmax180_);

  log << Log::INFO << "Sum of histJetAzimuth_pTmax75_100 bin heights after normalization: "
      << histJetAzimuth_pTmax75_100->sumBinHeights() << endl;
  log << Log::INFO << "Sum of histJetAzimuth_pTmax100_130 bin heights after normalization: "
      << histJetAzimuth_pTmax100_130->sumBinHeights() << endl;
  log << Log::INFO << "Sum of histJetAzimuth_pTmax130_180 bin heights after normalization: "
      << histJetAzimuth_pTmax130_180->sumBinHeights() << endl;
  log << Log::INFO << "Sum of histJetAzimuth_pTmax180_ bin heights after normalization: "
      << histJetAzimuth_pTmax180_->sumBinHeights() << endl;
}
