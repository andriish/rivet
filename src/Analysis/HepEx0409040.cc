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
  // Use histogram auto-booking
  histJetAzimuth_pTmax75_100  = bookHistogram1D(1, 2, 1, "Jet Jet azimuthal angle, pTmax=75..100");
  histJetAzimuth_pTmax100_130 = bookHistogram1D(2, 2, 1, "Jet Jet azimuthal angle, pTmax=100..130");
  histJetAzimuth_pTmax130_180 = bookHistogram1D(3, 2, 1, "Jet Jet azimuthal angle, pTmax=130..180");
  histJetAzimuth_pTmax180_    = bookHistogram1D(4, 2, 1, "Jet Jet azimuthal angle, pTmax>180");
}


// Do the analysis
void HepEx0409040::analyze(const Event & event) {
  Log& log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;

  // Analyse and print some info  
  const D0RunIIconeJets& jetpro = event.applyProjection(conejets);
  log << Log::INFO << "Jet multiplicity before any pT cut = " << jetpro.getNJets() << endl;
   
  // Find vertex and check  that its z-component is < 50 cm from the nominal IP
  //const PVertex& PV = event.applyProjection(p_vertex);
  /// @todo SEGV: either the HepMC event record is not filled properly or the F77-Wrapper functions are faulty
  /// @todo z- value assumed to be in mm, PYTHIA convention: dangerous!
  //if (fabs(PV().position().z())< 500.0) {
    list<HepEntity>::iterator jetpTmax = jetpro.jets->end();
    list<HepEntity>::iterator jet2ndpTmax = jetpro.jets->end();
    log << Log::DEBUG << "jetlist size = " << jetpro.jets->size() << endl;
    
    int Njet=0;
    for (list<HepEntity>::iterator jt = jetpro.jets->begin(); jt != jetpro.jets->end(); ++jt) {
      log << Log::DEBUG << "List item pT = " << jt->pT() << " E=" << jt->E << " pz=" << jt->pz << endl;
      if (jt->pT() > 40.0) ++Njet;
      log << Log::DEBUG << "Jet pT =" << jt->pT() << " y=" << jt->y() << " phi=" << jt->phi() << endl; 
      if (jetpTmax == jetpro.jets->end() || jt->pT() > jetpTmax->pT()) {
        jet2ndpTmax = jetpTmax;
        jetpTmax = jt;
      } else if (jet2ndpTmax == jetpro.jets->end() || jt->pT() > jet2ndpTmax->pT()) {
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
        log << Log::DEBUG << "Jet eta and pT requirements fulfilled" << endl;
        /// @todo Declare this eta cut via Analysis::addCut()?
        double etaMax = 3.0; //D0 calorimeter boundary
        bool addMuons = false; //Muons pass calorimeter almost without energy loss
        p_calmet.initialize(etaMax, addMuons);
        const CalMET& CaloMissEt = event.applyProjection(p_calmet);
        log << Log::DEBUG << "CaloMissEt.MET() = " << CaloMissEt.MET() << endl;
        if (CaloMissEt.MET() < 0.7*jetpTmax->pT()) {
	  
          double dphi = delta_phi(jetpTmax->phi(), jet2ndpTmax->phi());
          
          if (jetpTmax->pT() > 75.0 && jetpTmax->pT() <= 100.0)
            histJetAzimuth_pTmax75_100->fill(dphi, event.weight());
          else if (jetpTmax->pT() > 100.0 && jetpTmax->pT() <= 130.0)
            histJetAzimuth_pTmax100_130->fill(dphi, event.weight());
          else if (jetpTmax->pT() > 130.0 && jetpTmax->pT() <= 180.0)
            histJetAzimuth_pTmax130_180->fill(dphi, event.weight());
          else if (jetpTmax->pT() > 180.0)
            histJetAzimuth_pTmax180_->fill(dphi, event.weight());
          
        } //CalMET
      } //jets N, pT
    } //jets y (raqpidity) 
    
    //} //z-vertex
  
  
  // Finished
  log << Log::DEBUG << "Finished analyzing" << endl;
}


// Finalize
void HepEx0409040::finalize() { 
  Log& log = getLog();

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
