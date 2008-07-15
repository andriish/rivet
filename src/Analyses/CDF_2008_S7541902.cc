// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S7541902.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"


namespace Rivet {


  void CDF_2008_S7541902::init() {
    /// @todo Use histogram auto-booking when the paper is in HepData.
    /*  for (int i = 0 ; i < 3 ; i++) {
      stringstream title;
      title << "Et of jet " << i + 1;
      _histJetEt[i]    = bookHistogram1D(i+1,1,1,title.str().c_str());
    }
    _histJetMult   = bookHistogram1D(4,1,1, "Inclusive Jet Multiplicity");
    */
    _histJetMult   = bookHistogram1D("histJetMult", "Inclusive Jet Multiplicity", 5,-0.5,4.5);
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



  void CDF_2008_S7541902::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;
   
    // Get final state particles in event  
    const FinalState& part = applyProjection<FinalState>(event, "FS");
    const ParticleVector& particles =  part.particles();

    // Make kinematic cuts
    FourMomentum electronP, neutrinoP;
    for (ParticleVector::const_iterator p = particles.begin(); p != particles.end(); ++p) {
      // Pick out final state electrons and electron neutrinos
      const unsigned int abs_id = abs(p->getPdgId());

      if (abs_id == ELECTRON || abs_id == NU_E) {
        bool fromW = false;
        GenVertex *startVtx = (p->getHepMCParticle()).production_vertex();       
        if (startVtx) {
          for (GenVertex::particle_iterator pIt = startVtx->particles_begin(HepMC::ancestors);
               pIt != startVtx->particles_end(HepMC::ancestors); ++pIt) {
            if (fabs((*pIt)->pdg_id())==24) fromW = true;
          }
        } else {
          log << Log::WARN << "No vertex found for final state particles" << endl;
        }
        if (fromW && abs_id == ELECTRON) electronP = p->getMomentum();
        if (fromW && abs_id == NU_E)     neutrinoP = p->getMomentum();
      }
    }

    /// @todo Declare these cuts
    bool fail = false;
    if (electronP.pT() < 20.0 || fabs(electronP.pseudorapidity()) > 1.1) fail = true;
    if (neutrinoP.pT() < 30.0)  fail = true;
    if (electronP.pT() == 0.)  log << Log::TRACE << "No final state electron found, probably outside eta range" << endl;
    if (neutrinoP.pT() == 0.)  log << Log::TRACE << "No final state neutrino found, probably outside eta range" << endl;
    
    if (!fail) {
      /// @todo Should add this as a function in FourMom
      float MT = electronP.pT()*neutrinoP.pT() - 
        electronP.px()*neutrinoP.px() -
        electronP.py()*neutrinoP.py();
      MT = sqrt(2.0*MT);
      if (MT < 20.) fail = true;
    }
    if (fail) vetoEvent(event);

    // JetClu 0.4 radius, merging fraction = 0.75 
    // Remove electron (and neutrinos/muons) before clustering
    const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
    const PseudoJets& jets = jetpro.getPseudoJetsByPt();
    vector<float> jetPt;  
    size_t njets = 0; 
    for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      if (fabs(jt->pseudorapidity()) < 2.0 ){
        if (jt->perp() > 20.) jetPt.push_back(jt->perp());
        if (jt->perp() > 25.) ++njets;
      } 
    } 
    for (size_t i = 0; i <= njets ; ++i) {
      _histJetMult->fill(i,event.weight());
    }
    for (size_t i = 0; i < 5; ++i){
      if (jetPt.size() > i) {
        _histJetEt[i]->fill(jetPt[i],event.weight());
      }
    }

  }
  
  
  
  // Finalize
  void CDF_2008_S7541902::finalize() { 
    float xsec = crossSection()/sumOfWeights();
    for (int i=0; i<4; i++) {
      _histJetEt[i]->scale(xsec);
    }
   _histJetMult->scale(xsec);
  }
  
    
}
