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
  _histJetAzimuth_pTmax75_100  = bookHistogram1D(1, 2, 1, "Jet Jet azimuthal angle, pTmax=75..100");
  _histJetAzimuth_pTmax100_130 = bookHistogram1D(2, 2, 1, "Jet Jet azimuthal angle, pTmax=100..130");
  _histJetAzimuth_pTmax130_180 = bookHistogram1D(3, 2, 1, "Jet Jet azimuthal angle, pTmax=130..180");
  _histJetAzimuth_pTmax180_    = bookHistogram1D(4, 2, 1, "Jet Jet azimuthal angle, pTmax>180");
}


// Do the analysis
void HepEx0409040::analyze(const Event & event) {
  Log& log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;

  // Analyse and print some info  
  const D0ILConeJets& jetpro = event.applyProjection(_conejetsproj);
  log << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.getNJets() << endl;
   
  // Find vertex and check  that its z-component is < 50 cm from the nominal IP
  //const PVertex& pv = event.applyProjection(_vertexproj);
  /// @todo SEGV: either the HepMC event record is not filled properly or the F77-Wrapper functions are faulty
  /// @todo z- value assumed to be in mm, PYTHIA convention: dangerous!
  //if (fabs(pv.getPrimaryVertex().position().z()) < 500.0) {
    //list<HepEntity>::const_iterator jetpTmax = jetpro._jets.end();
    //list<HepEntity>* jets = (const list<HepEntity>*) 
    const list<HepEntity>* jets = jetpro.getJets();
    list<HepEntity>::const_iterator jetpTmax = jets->end();
    list<HepEntity>::const_iterator jet2ndpTmax = jets->end();
    log << Log::DEBUG << "jetlist size = " << jets->size() << endl;

    int Njet = 0;
    for (list<HepEntity>::const_iterator jt = jets->begin(); jt != jets->end(); ++jt) {
      log << Log::DEBUG << "List item pT = " << jt->pT() << " E=" << jt->E << " pz=" << jt->pz << endl;
      if (jt->pT() > 40.0) ++Njet;
      log << Log::DEBUG << "Jet pT =" << jt->pT() << " y=" << jt->y() << " phi=" << jt->phi() << endl; 
      if (jetpTmax == jets->end() || jt->pT() > jetpTmax->pT()) {
        jet2ndpTmax = jetpTmax;
        jetpTmax = jt;
      } else if (jet2ndpTmax == jets->end() || jt->pT() > jet2ndpTmax->pT()) {
        jet2ndpTmax = jt;
      }
    }
    
    if (Njet >= 2) {
      log << Log::DEBUG << "Jet multiplicity after pT > 40GeV cut = " << Njet << endl; 
    }

    if (jets->size()>=2 && jet2ndpTmax->pT() > 40.) {
      if (fabs(jetpTmax->y())<0.5 && fabs(jet2ndpTmax->y())<0.5) {
        log << Log::DEBUG << "Jet eta and pT requirements fulfilled" << endl;
        /// @todo Should this commented eta cut be happening via a FinalState configuration?
        /// @todo Declare this eta cut via Analysis::addCut()?
        //double etaMax = 3.0; //D0 calorimeter boundary
        const TotalVisibleMomentum& caloMissEt = event.applyProjection(*_calmetproj);
        log << Log::DEBUG << "CaloMissEt.getMomentum().perp() = " << caloMissEt.getMomentum().perp() << endl;
        if (caloMissEt.getMomentum().perp() < 0.7*jetpTmax->pT()) {
          double dphi = delta_phi(jetpTmax->phi(), jet2ndpTmax->phi());
          
          if (jetpTmax->pT() > 75.0 && jetpTmax->pT() <= 100.0)
            _histJetAzimuth_pTmax75_100->fill(dphi, event.weight());
          else if (jetpTmax->pT() > 100.0 && jetpTmax->pT() <= 130.0)
            _histJetAzimuth_pTmax100_130->fill(dphi, event.weight());
          else if (jetpTmax->pT() > 130.0 && jetpTmax->pT() <= 180.0)
            _histJetAzimuth_pTmax130_180->fill(dphi, event.weight());
          else if (jetpTmax->pT() > 180.0)
            _histJetAzimuth_pTmax180_->fill(dphi, event.weight());
          
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
  normalize(_histJetAzimuth_pTmax75_100);
  normalize(_histJetAzimuth_pTmax100_130);
  normalize(_histJetAzimuth_pTmax130_180);
  normalize(_histJetAzimuth_pTmax180_);

  log << Log::INFO << "Sum of histJetAzimuth_pTmax75_100 bin heights after normalization: "
      << _histJetAzimuth_pTmax75_100->sumBinHeights() << endl;
  log << Log::INFO << "Sum of histJetAzimuth_pTmax100_130 bin heights after normalization: "
      << _histJetAzimuth_pTmax100_130->sumBinHeights() << endl;
  log << Log::INFO << "Sum of histJetAzimuth_pTmax130_180 bin heights after normalization: "
      << _histJetAzimuth_pTmax130_180->sumBinHeights() << endl;
  log << Log::INFO << "Sum of histJetAzimuth_pTmax180_ bin heights after normalization: "
      << _histJetAzimuth_pTmax180_->sumBinHeights() << endl;
}
