// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S7541902.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"


namespace Rivet {


  // Book histograms
  void CDF_2008_S7541902::init() {
    /// @todo Use histogram auto-booking when the paper is in HepData.
    /// @todo get root output
    _histJetMult   = bookHistogram1D("histJetMult", "Inclusive Jet Multiplicity", 5,-0.5,4.5);
    // _histJetEt[0]    = bookHistogram1D(1,1,1,"First Jet Et");
    const int bins1 = 13;
    const double xbins1[bins1] = {20.,25.,30.,35.,40.,50.,60.,75.,90.,110.,150.,195.,350.};
    _histJetEt[0]    = bookHistogram1D("histJetEt_1jet", "First Jet Et", vector<double>(xbins1, xbins1+bins1));
    const int bins2 = 10;
    const double xbins2[bins2] = {20.,25.,30.,35.,40.,50.,60.,75.,95.,190.};
    _histJetEt[1]  = bookHistogram1D("histJetEt_2jet", "Second Jet Et", vector<double>(xbins2, xbins2+bins2));
    const int bins3 = 6;
    const double xbins3[bins3] = {20.,25.,30.,35.,45.,80.};
    _histJetEt[2]  = bookHistogram1D("histJetEt_3jet", "Third Jet Et", vector<double>(xbins3, xbins3+bins3));
    const int bins4 = 3;
    const double xbins4[bins4] = {20.,25.,35.};
    _histJetEt[3]  = bookHistogram1D("histJetEt_4jet", "Fourth Jet Et", vector<double>(xbins4, xbins4+bins4));
  }

  // Do the analysis
  void CDF_2008_S7541902::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;
    
    bool fail = false;
    // Get final state particles in event  
    const FinalState& part = event.applyProjection(_fsproj);
    const ParticleVector& particles =  part.particles();

    // Make kinematic cuts
    FourMomentum electronP(0.,0.,0.,0.);
    FourMomentum neutrinoP(0.,0.,0.,0.); 
    for (ParticleVector::const_iterator p = particles.begin(); p != particles.end(); ++p) {
      // pick out final state electrons and electron neutrinos
      int abs_id = abs(p->getPdgId());
      if (abs_id == 11 || abs_id == 12) {
        bool fromW = false;
	GenVertex *startVtx=(p->getHepMCParticle()).production_vertex();       
	if(startVtx!=0) {
	  for(GenVertex::particle_iterator pIt = startVtx->particles_begin(HepMC::ancestors);
	      pIt != startVtx->particles_end(HepMC::ancestors); ++pIt){
	    if(fabs((*pIt)->pdg_id())==24) fromW = true;
	  }
        
	} else log << Log::WARN << "No vertex found for final state particles" << endl;
        if (fromW && abs_id == 11)  electronP = p->getMomentum();
        if (fromW && abs_id == 12)  neutrinoP = p->getMomentum();
      }
    }
    if (electronP.pT() < 20. || fabs(electronP.pseudorapidity()) > 1.1) fail = true;
    if (neutrinoP.pT() < 30.)  fail = true;
    if (electronP.pT() == 0.)  log << Log::WARN << "No final state electron found" << endl;
    if (neutrinoP.pT() == 0.)  log << Log::WARN << "No final state neutrino found" << endl;
    if (!fail) {
      // should add this as a function in FourMom
      float MT = electronP.pT()*neutrinoP.pT() - 
	electronP.px()*neutrinoP.px() -
	electronP.py()*neutrinoP.py();
      MT = sqrt(2.0*MT);
      if(MT < 20.) fail = true;
    }
    // if the event passes the kinematic cuts
    if (!fail) {      
      // JetClu 0.4 radius, merging fraction = 0.75 
      // Remove electron (and neutrinos/muons) before clustering
      const FastJets& jetpro = event.applyProjection(_conejetsproj);
      // get pT ordered pseudo jets (whats that?)
      const PseudoJets& jets = jetpro.getPseudoJetsPt();
      vector<float> jetPt;  
      int njets = 0;      
      for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
 	if (fabs(jt->pseudorapidity()) < 2.0 ){
	  if(jt->perp() > 20.) jetPt.push_back(jt->perp());
	  if(jt->perp() > 25.)  ++njets;
	} 
      } 
      for (int i=0; i<=njets ; i++) {
 	_histJetMult->fill(i,event.weight());
      }
      for (unsigned int i=0; i<5; i++){
      	if(jetPt.size() > i) _histJetEt[i]->fill(jetPt[i],event.weight());
      }
    } // end passes
    // Finished
    log << Log::DEBUG << "Finished analyzing" << endl;
    
  }
  
  
  
  // Finalize
  void CDF_2008_S7541902::finalize() { 
    float xsec = crossSection()/sumOfWeights();
    _histJetMult->scale(xsec);
    for (int i=0; i<4; i++) _histJetEt[i]->scale(xsec);
  }
  
  
}
