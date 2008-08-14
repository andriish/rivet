// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S7541902.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/AnalysisHandler.hh"
#include <algorithm>

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

    //get the hardest electron in the event
    
    const ChargedLeptons &chLeptons = applyProjection<ChargedLeptons>(event, "ChLeptons");
    ParticleVector theLeptons = chLeptons.chargedLeptons();
    sort(theLeptons.begin(), theLeptons.end(), Particle::byETDescending());
    
    Particle electron;
    bool gotLepton = false;
    bool stop = false;
    for(ParticleVector::const_iterator p = theLeptons.begin();
        !stop && !gotLepton && p != theLeptons.end(); ++p){
      if(p->momentum().Et() > _electronETCut){
        if(abs(p->getPdgId())==11){
          electron = *p;
          gotLepton = true;
        }
      }else{
        stop = true;
      }
    }

    if(!gotLepton){
      vetoEvent(event);
      return;
    }
    
    // get the neutrino
    
    FourMomentum neutrinoP;
    bool gotNeutrino=false;
    
    for(ParticleVector::const_iterator p = particles.begin(); 
        !gotNeutrino && p != particles.end(); ++p){
      const unsigned int abs_id = abs(p->getPdgId());
      
      if(abs_id == NU_E){
        neutrinoP = p->momentum();
        if(neutrinoP.Et() > _eTmissCut) gotNeutrino = true;
      }
    }
    
    if(!gotNeutrino){
      vetoEvent(event);
      return;
    }
    
    // get the jets
    
    const FastJets &jetProj = applyProjection<FastJets>(event, "Jets");
    
    Jets theJets = jetProj.getJetsByPt(15.0);
    sort(theJets.begin(), theJets.end(), Particle::byETDescending());
    
    Jets foundJets;
    stop = false;
    
    for(Jets::const_iterator jIt =theJets.begin(); 
        !stop && foundJets.size()!=4 && jIt != theJets.end(); ++jIt){
      
      if(jIt->momentum().Et()>_jetEtCut){
        if(! jIt->containsParticle(electron))
          foundJets.push_back(*jIt);
      }else{
        stop = true;
      }
    }
    
    double mT2 = electron.momentum().pT() * neutrinoP.pT() - 
    electron.momentum().px() * neutrinoP.px() -
    electron.momentum().py() * neutrinoP.py();
           
    if(mT2 < _mT2Cut ){
      vetoEvent(event);
      return;
    }
    
    for(size_t i = 0; i != foundJets.size(); ++i){
      _histJetEt[i]->fill(foundJets[i].momentum().Et(),event.weight());
    }
    
    return;
  }
  
  
  
  // Finalize
  void CDF_2008_S7541902::finalize() { 
    float xsec = crossSection()/getHandler().sumOfWeights();
    for (int i=0; i<4; i++) {
      _histJetEt[i]->scale(xsec);
    }
   _histJetMult->scale(xsec);
  }
  
    
}
