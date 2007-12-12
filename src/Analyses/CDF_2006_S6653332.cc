// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2006_S6653332.hh"
#include "Rivet/RivetAIDA.hh"
#include "HepPDT/ParticleID.hh"
//#include "Rivet/Tools/Logging.hh"


namespace Rivet {


  /// Analysis dependent cuts on vertex tracks in SVertex projection 
  /// Since the analysis specific cuts are very complex, they are not 
  /// implemented in the projection and are instead passed via a function (object).
  /// SVertex member function implementation below 
  /// in: reference to instance of SVertex projection, ParticleVector of
  ///     vertex to be analyzed, primary (Gen)Vertex
  /// out: FourMomentum = visible Momentum of vertex (selected tracks), 
  /// return bool: cuts passed? 1 : 0 

  bool CDF_2006_S6653332::applyVtxTrackCuts(SVertex& svtx, ParticleVector& chfsvtx, 
				       const HepMC::GenVertex& gpvtx,
				       FourMomentum vtxVisMom) {
    ///Check vertex final state charged particles, if fulfilling track criteria
    int pass1trk1pTdcaSig25 = 0;
    int pass1trk05pTdcaSig25 = 0;
    int pass2trk15pTdcaSig3 = 0;
    int pass2trk1pTdcaSig3 = 0;
 
    //Log log = getLog();

    for (size_t i=0; i<chfsvtx.size(); ++i) {
      double IPsig = svtx.get2dDCAsig(chfsvtx[i].getHepMCParticle(), gpvtx);

      //log << Log::DEBUG << "CDF_2006_S6653332::applyVtxTrackCuts: IPsig = " 
      //  << IPsig << endl;

      if (chfsvtx[i].getMomentum().pT() > 0.5) 
        vtxVisMom += chfsvtx[i].getMomentum();

      if (chfsvtx.size()>=3 && IPsig>2.5) { // 1st pass
        if (chfsvtx[i].getMomentum().pT() > 1.0) pass1trk1pTdcaSig25++;
        else if (chfsvtx[i].getMomentum().pT()>0.5) pass1trk05pTdcaSig25++;
      }
      if (chfsvtx.size()>=2 && IPsig>3.) { // 2nd pass
        if (chfsvtx[i].getMomentum().pT() > 1.5) pass2trk15pTdcaSig3++;
        else if (chfsvtx[i].getMomentum().pT()>1.) pass2trk1pTdcaSig3++;
      } 
    }

    // log << Log::DEBUG << "CDF_2006_S6653332::applyVtxTrackCuts: pass1trk1pTdcaSig25=" 
    // << pass1trk1pTdcaSig25 << " pass1trk05pTdcaSig25="  
    // << pass1trk05pTdcaSig25 << " pass2trk15pTdcaSig3=" << pass2trk15pTdcaSig3 
    // << " pass2trk1pTdcaSig3=" << pass2trk1pTdcaSig3 << endl;     
    
    bool cutspassed = false;
    if (pass1trk1pTdcaSig25>=1 && pass1trk1pTdcaSig25+pass1trk05pTdcaSig25>=3 ||
        pass2trk15pTdcaSig3>=1 && pass2trk15pTdcaSig3+pass2trk1pTdcaSig3>=2) {
      cutspassed = true;
      //log << Log::DEBUG << "CDF_2006_S6653332::applyVtxTrackCuts: Vertex track cuts passed!" << endl;
    }
    
    return cutspassed;
  }
  


  // Book histograms
  void CDF_2006_S6653332::init() {
    /// @todo Use histogram auto-booking if there are comparison datasets in HepData.
    // Will be replaced as soon as HepData is updated and rivet/data subdirectory
    // contains relevant data histo's
    _histJetsPt    = bookHistogram1D("JetsPt", "All jets Pt", 10, 0., 100.);
    _histJetsEta   = bookHistogram1D("JetsEta", "All jets Pseudorapidity", 20, -2., 2.);
    _histbJetsPt   = bookHistogram1D("bJetsPt", "b jets Pt", 10, 0., 100.);
    _histbJetsEta  = bookHistogram1D("bJetsEta", "b jets Pseudorapidity", 20, -2., 2.);
  }



  // Do the analysis
  void CDF_2006_S6653332::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;
    
    // Find primary vertex 
    const PVertex& pv = event.applyProjection(_pvtxproj);
    if (fabs(pv.getPrimaryVertex().position().z()) < 600.0) {

    const TotalVisibleMomentum& caloMissEt = event.applyProjection(_calmetproj);
    log << Log::DEBUG << "CaloMissEt.getMomentum().pT() = " << caloMissEt.getMomentum().pT() << endl;
    
    if (caloMissEt.getMomentum().pT()<25. && caloMissEt.getSET()<150.) {
      
      // Get leptons in event
      const ChargedLeptons& chlep = event.applyProjection(_chleproj);
      const ParticleVector& leps = chlep.chargedLeptons();
      
      // Sort 2 leading leptons in pT
      vector<ParticleVector::const_iterator> dilep;
      for (ParticleVector::const_iterator p = leps.begin(); p != leps.end(); ++p) {
        if (dilep.size()==0) dilep.push_back(p);
        else if (dilep.size()==1) {
          if (p->getMomentum().pT() > dilep[0]->getMomentum().pT() )  { 
            dilep.push_back(dilep[0]); 
            dilep[0]= p; 
          }
          else dilep.push_back(p);
        }
        else { // dilep.size == 2
          if (p->getMomentum().pT() > dilep[0]->getMomentum().pT() )  { 
            dilep[1] = dilep[0];
            dilep[0] = p;
          }
          else if (p->getMomentum().pT() > dilep[1]->getMomentum().pT() )  { 
            dilep[1] = p; 
          }
        }
      }
      
      
      // Require exactly 2 leptons, same flavour, opposite charge, pTmin > 10 GeV.
      if (dilep.size()==2 && dilep[1]->getMomentum().pT() > 10. 
          && dilep[0]->getPdgId() == -dilep[1]->getPdgId() ) {

        //cout << "CDF_2006_S6653332::analyze: Z pre-cuts passed!" << endl;
        
        // Muon and electron geometrical acceptance requirement
        /// @todo Use ParticleName enum for clarity.
        if ((abs(dilep[0]->getPdgId() == 11) && fabs(dilep[0]->getMomentum().pseudorapidity()) < 3.6 && fabs(dilep[1]->getMomentum().pseudorapidity()) < 3.6) ||
            (abs(dilep[0]->getPdgId() == 13) && fabs(dilep[0]->getMomentum().pseudorapidity()) < 1.5 && fabs(dilep[1]->getMomentum().pseudorapidity()) < 1.5)    
            ) {
          
          // Dilepton invariant mass cut
          FourMomentum ll = dilep[0]->getMomentum() + dilep[1]->getMomentum();
          log << Log::DEBUG << "ll.mass()=" << ll.mass() << endl;
          /// @todo Use const doubles insted of all these magic numbers...
          if (66. < ll.mass()  &&  ll.mass() < 116.) {
	    
            log << Log::DEBUG << "Invariant di-lepton mass ll.mass() = " << ll.mass() << endl;
            
            double rapmax;
            /// @todo Use ParticleName enum
            if (abs(dilep[0]->getPdgId()) == 11) rapmax = 1.1; //electron
            else rapmax = 1.0; //muon
            
            // More restrictive pT requirement on (one) trigger lepton
            if ( (dilep[0]->getMomentum().pT() > 18.0 && 
                  fabs(dilep[0]->getMomentum().pseudorapidity()) < rapmax) || 
                 (dilep[1]->getMomentum().pT() > 18.0 && 
                  fabs(dilep[1]->getMomentum().pseudorapidity()) < rapmax)) {
              
              
              const VetoedFinalState& vfs = event.applyProjection(_vfsproj);
              // Check if leptons are isolated
              double R_isol=0.4;
              vector<double> Ehad(2, 0.); //dim=2 for 2 leptons
              vector<double> Eem(2, 0.); //dim=2 for 2 leptons
              double fh = 0.8; // Hadronic calo. fraction of charged hadrons (simplified assumption)
              bool l_isol = true; // Both leptons have to be isolated 
              // => one veto variable suffices
              for (ParticleVector::const_iterator p = vfs.particles().begin(); 
                   p != vfs.particles().end(); ++p) {
                
                HepPDT::ParticleID pInfo = p->getPdgId();
                
                // Determine energy depositions around lepton(s)
                /// @todo Use ParticleName enum.
                for (int lind=0; lind<2; ++lind) {
                  if (deltaR(p->getMomentum(), dilep[lind]->getMomentum()) < R_isol) {
                    if (abs(p->getPdgId()) == 11 
                        || abs(p->getPdgId()) == 13 
                        || abs(p->getPdgId()) == 22) {
                      Eem[lind] += p->getMomentum().E();
                    } else if (pInfo.threeCharge() != 0) {//charged hadron: had + em calo fraction
                      Eem[lind] += (1.-fh)*p->getMomentum().E();
                      Ehad[lind] += fh*p->getMomentum().E();
                    } else { // Neutral hadron assumed to have vanishing EM fraction
                      Ehad[lind] += p->getMomentum().E();
                    }
                  }
                }
              }
              
              /// @todo This is cryptic: use some constants instead of all these magic numbers.
              /// @todo Do leptons fulfill isolation criteria?
              for (int lind=0; lind < 2; ++lind) {
                if (abs(dilep[lind]->getPdgId())==11) { // electron
                  if (!(Eem[lind]>0.)) l_isol = false;
                  else if (Ehad[lind]/Eem[lind] > 0.055+0.00045*dilep[lind]->getMomentum().E()) l_isol = false;
                }
                else { //muon
                  if (dilep[lind]->getMomentum().E() <= 100.) {
                    if (Ehad[lind]-dilep[lind]->getMomentum().E() > 6.
                        || Eem[lind]-dilep[lind]->getMomentum().E() > 2.) l_isol = false;
                  }
                  else { // p > 100.
                    if (Ehad[lind]-dilep[lind]->getMomentum().E() > // E_had = 0.0280 
                        6.+(dilep[lind]->getMomentum().E()-100.)*0.0280
                        || Eem[lind]-dilep[lind]->getMomentum().E() > // E_em = 0.0115 
                        2.+(dilep[lind]->getMomentum().E()-100.)*0.0115) l_isol = false;
                  }
                }
              }
              
              if (l_isol) { // Both leptons are isolated
                
                log << Log::DEBUG << "Both leptons are isolated" << endl;
                
                const FastJets& jetpro = event.applyProjection(_jetsproj);
                /// @todo Don't expose FastJet objects in Rivet analyses.
                const PseudoJets& jets = jetpro.getPseudoJetsPt();
                
                // Remove leptons from list of jets
                _jetaxes.clear();
                for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
                  const FourMomentum j(jt->E(), jt->px(), jt->py(), jt->pz());
                  if (deltaR(dilep[0]->getMomentum(), j) > R_isol && 
                      deltaR(dilep[1]->getMomentum(), j) > R_isol) { _jetaxes.push_back(j); }
                }
                
                
                // Fill jets Et and eta distributions.
                for (size_t i=0; i < _jetaxes.size(); ++i) {
                  if (_jetaxes[i].Et() > 20. && 
                      fabs(_jetaxes[i].pseudorapidity()) < 1.5) {
                    _histJetsPt->fill(_jetaxes[i].Et(), event.weight());
                    _histJetsEta->fill(_jetaxes[i].pseudorapidity(), event.weight());
                  }
                }
                
                // Determine secondary vertices
                const SVertex& svtx =  event.applyProjection(_svtxproj);
                
                const vector<FourMomentum>& taggedJets = svtx.getTaggedJets();
                
                log << Log::DEBUG << "taggedJets.size()=" << taggedJets.size() << endl;

                // Fill tagged jets Et and eta distributions
                for (size_t i=0; i < taggedJets.size(); ++i) {
                  if (taggedJets[i].Et() > 20. && 
                      fabs(taggedJets[i].pseudorapidity()) < 1.5) {
                    _histJetsPt->fill(taggedJets[i].Et(), event.weight());
                    _histJetsEta->fill(taggedJets[i].pseudorapidity(), event.weight());
                  }
                }
                
              } //both leptons are isolated
            } //dilepton invariant mass cut
          } //mu,e geometrical acceptance requirement
        } //at least one trigger lepton (eta, pT) requirement	
      } //opposite charged leptons, pT>10
      
    }//Cal. missing Et and scalar Et (Ht) cut
    
    }// primary vertex z-position requirement
    
    // Finished
    log << Log::DEBUG << "Finished analyzing" << endl;
  }
  

  
  // Finalize
  void CDF_2006_S6653332::finalize() { 
    normalize(_histJetsPt);
    normalize(_histJetsEta);
    normalize(_histbJetsPt);
    normalize(_histbJetsEta);
  }


}
