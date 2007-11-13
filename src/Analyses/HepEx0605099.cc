// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/HepEx0605099.hh"
#include "Rivet/RivetAIDA.hh"
#include "HepPDT/ParticleID.hh"
//#include "Rivet/Tools/Logging.hh"



namespace Rivet {



  /// Analysis dependend cuts on vertex tracks in SVertex Projection 
  /// Since very complex analysis specific cuts, not implemented in Projection
  /// source code to keep it universal. 
  /// SVertex member function implementation below 
  /// in: reference to instance of SVertex projection, ParticleVector of
  ///     vertex to be analyzed, primary (Gen)Vertex
  /// out: LorentzVector = visible Momentum of vertex (selected tracks), 
  /// return bool: cuts passed? 1 : 0 

  bool HepEx0605099::applyVtxTrackCuts(SVertex& svtx, ParticleVector& chfsvtx, 
				       const HepMC::GenVertex& gpvtx,
				       LorentzVector vtxVisMom) {
    ///Check vertex final state charged particles, if fulfilling track criteria
    int pass1trk1pTdcaSig25 = 0;
    int pass1trk05pTdcaSig25 = 0;
    int pass2trk15pTdcaSig3 = 0;
    int pass2trk1pTdcaSig3 = 0;
 
    //Log log = getLog();

    for (unsigned int i=0; i<chfsvtx.size(); ++i) {
      double IPsig = svtx.get2dDCAsig(chfsvtx[i].getHepMCParticle(), gpvtx);

      //log << Log::DEBUG << "HepEx0605099::applyVtxTrackCuts: IPsig = " 
      //  << IPsig << endl;

      if (chfsvtx[i].getMomentum().perp()>0.5) 
	vtxVisMom += chfsvtx[i].getMomentum();

      if (chfsvtx.size()>=3 && IPsig>2.5) {// 1st pass
	if (chfsvtx[i].getMomentum().perp()>1.) pass1trk1pTdcaSig25++;
	else if (chfsvtx[i].getMomentum().perp()>0.5) pass1trk05pTdcaSig25++;
      }
      if (chfsvtx.size()>=2 && IPsig>3.) {// 2nd pass
	if (chfsvtx[i].getMomentum().perp()>1.5) pass2trk15pTdcaSig3++;
	else if (chfsvtx[i].getMomentum().perp()>1.) pass2trk1pTdcaSig3++;
      }      
    }

    //log << Log::DEBUG << "HepEx0605099::applyVtxTrackCuts: pass1trk1pTdcaSig25=" 
    //<< pass1trk1pTdcaSig25 << " pass1trk05pTdcaSig25="  
    //<< pass1trk05pTdcaSig25 << " pass2trk15pTdcaSig3=" << pass2trk15pTdcaSig3 
    //<< " pass2trk1pTdcaSig3=" << pass2trk1pTdcaSig3 << endl;     
    
    bool cutspassed = false;
    if (pass1trk1pTdcaSig25>=1 && pass1trk1pTdcaSig25+pass1trk05pTdcaSig25>=3 ||
	pass2trk15pTdcaSig3>=1 && pass2trk15pTdcaSig3+pass2trk1pTdcaSig3>=2) {
      cutspassed = true;
      //log << Log::DEBUG << "HepEx0605099::applyVtxTrackCuts: Vertex track cuts passed!" << endl;
    }
    
    return cutspassed;
  }
  




  // Book histograms
  void HepEx0605099::init() {
    // Using histogram auto-booking is preferable if there are comparison datasets in HepData.
    // Will be replaced as soon as HepData is updated and rivet/data subdirectory
    // contains relevant data histo's
    _histJetsPt         = bookHistogram1D("JetsPt", "All jets Pt", 10, 0., 100.);
    _histJetsEta       = bookHistogram1D("JetsEta", "All jets Pseudorapidity", 20, -2., 2.);
    _histbJetsPt         = bookHistogram1D("bJetsPt", "b jets Pt", 10, 0., 100.);
    _histbJetsEta       = bookHistogram1D("bJetsEta", "b jets Pseudorapidity", 20, -2., 2.);

  }




  // Do the analysis
  void HepEx0605099::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;
    
    // Find primary vertex 
    const PVertex& pv = event.applyProjection(_pvtxproj);
    if (fabs(pv.getPrimaryVertex().position().z()) < 600.0) {

    const TotalVisibleMomentum& caloMissEt = event.applyProjection(_calmetproj);
    log << Log::DEBUG << "CaloMissEt.getMomentum().perp() = " << caloMissEt.getMomentum().perp() << endl;
    
    if (caloMissEt.getMomentum().perp()<25. && caloMissEt.getSET()<150.) {
      
      /// get leptons in event
      const ChargedLeptons& chlep = event.applyProjection(_chleproj);
      const ParticleVector& leps = chlep.chargedLeptons();

      /// sort 2 leading leptons in pT
      vector<ParticleVector::const_iterator> dilep;
      for (ParticleVector::const_iterator p = leps.begin(); p != leps.end(); ++p) {
	if (dilep.size()==0) dilep.push_back(p);
	else if (dilep.size()==1) {
	  if (p->getMomentum().perp() > dilep[0]->getMomentum().perp() )  { 
	    dilep.push_back(dilep[0]); 
	    dilep[0]= p; 
	  }
	  else dilep.push_back(p);
	}
	else { //dilep.size==2
	  if (p->getMomentum().perp() > dilep[0]->getMomentum().perp() )  { 
	    dilep[1] = dilep[0];
	    dilep[0] = p;
	  }
	  else if (p->getMomentum().perp() > dilep[1]->getMomentum().perp() )  { 
	    dilep[1] = p; 
	  }
	}
      }

      
      ///require exactly 2 leptons, same flavour, opposite charge, pTmin>10GeV
      if (dilep.size()==2 && dilep[1]->getMomentum().perp()>10. 
	  && dilep[0]->getPdgId() == -dilep[1]->getPdgId() ) {

	//cout << "HepEx0605099::analyze: Z pre-cuts passed!" << endl;

	/// Muon and electron geometrical acceptance requirement
	if ((abs(dilep[0]->getPdgId()==11) && fabs(dilep[0]->getMomentum().pseudoRapidity())<3.6 && fabs(dilep[1]->getMomentum().pseudoRapidity())<3.6) ||
	    (abs(dilep[0]->getPdgId()==13) && fabs(dilep[0]->getMomentum().pseudoRapidity())<1.5 && fabs(dilep[1]->getMomentum().pseudoRapidity())<1.5)    
	    ) {
	  
	  ///dilepton invariant mass cut
	  LorentzVector ll = dilep[0]->getMomentum() + dilep[1]->getMomentum();
	  cout << "ll.m()=" << ll.m() << endl;
	  if (66. < ll.m()  &&  ll.m() < 116.) {
	    
	    log << Log::DEBUG << "Invariant di-lepton mass ll.m() = " << ll.m() << endl;
	    
	    double rapmax;
	    if (abs(dilep[0]->getPdgId())==11) rapmax = 1.1; //electron
	    else rapmax = 1.0; //muon
	    
	    ///more restrictive pT requirement on (one) trigger lepton
	    if ( (dilep[0]->getMomentum().perp()>18.
		  && fabs(dilep[0]->getMomentum().pseudoRapidity()) < rapmax)
		 || (dilep[1]->getMomentum().perp()>18. 
		     && fabs(dilep[1]->getMomentum().pseudoRapidity()) < rapmax)) {
	      
	      
	      const VetoedFinalState& vfs = event.applyProjection(_vfsproj);
	      ///check if leptons are isolated
	      double R_isol=0.4;
	      vector<double> Ehad(2, 0.); //dim=2 for 2 leptons
	      vector<double> Eem(2, 0.); //dim=2 for 2 leptons
	      double fh = 0.8; //hadronic calo. fraction of charged hadrons (simplified assumption)
	      bool l_isol = true; //both leptons have to be isolated 
	      // => one veto variable suffices
	      for (ParticleVector::const_iterator p = vfs.particles().begin(); 
		   p != vfs.particles().end(); ++p) {
		
		HepPDT::ParticleID pInfo = p->getPdgId();
		
		///determine energy depositions around lepton(s)
		for (int lind=0; lind<2; ++lind) {
		  if (p->getMomentum().deltaR(dilep[lind]->getMomentum()) < R_isol) {
		    if (abs(p->getPdgId())==11 
			|| abs(p->getPdgId())==13 
			|| abs(p->getPdgId())==22) Eem[lind] += p->getMomentum().e();
		    //else if (p->getPdgId().threeCharge() != 0) {
		    else if (pInfo.threeCharge() != 0) {//charged hadron: had + em calo fraction
		      Eem[lind] += (1.-fh)*p->getMomentum().e();
		      Ehad[lind] += fh*p->getMomentum().e();
		    }
		    else //neutral hadron assumed to have vanishing EM fraction
		      Ehad[lind] += p->getMomentum().e();
		  }
		}
	      }
	      
	      ///do leptons fulfill isolation criteria?
	      for (int lind=0; lind<2; ++lind) {
		if (abs(dilep[lind]->getPdgId())==11) { //electron
		  if (!(Eem[lind]>0.)) l_isol = false;
		  else if (Ehad[lind]/Eem[lind] > 0.055+0.00045*dilep[lind]->getMomentum().e()) l_isol = false;
		}
		else { //muon
		  if (dilep[lind]->getMomentum().e() <= 100.) {
		    if (Ehad[lind]-dilep[lind]->getMomentum().e() > 6.
			|| Eem[lind]-dilep[lind]->getMomentum().e() > 2.) l_isol = false;
		  }
		  else { // p > 100.
		    if (Ehad[lind]-dilep[lind]->getMomentum().e() > //E_had=0.0280    
			6.+(dilep[lind]->getMomentum().e()-100.)*0.0280
			|| Eem[lind]-dilep[lind]->getMomentum().e() > //E_em=0.0115    
			2.+(dilep[lind]->getMomentum().e()-100.)*0.0115) l_isol = false;
		  }
		}
	      }
	      
	      if (l_isol) { //both leptons are isolated

		log << Log::DEBUG << "Both leptons are isolated" << endl;


		/// @todo Remove compile-time flags like this when possible.
                #ifdef HAVE_FASTJET
	        const FastJets& jetpro = event.applyProjection(_jetsproj);
                #else
	        const D0ILConeJets& jetpro = event.applyProjection(_jetsproj);
                #endif
	      
                /// @todo Remove compile-time flags like this when possible.
                #ifdef HAVE_FASTJET
                typedef vector<fastjet::PseudoJet> Jets;
                const Jets& jets = jetpro.getJetsPt();
                #else
                typedef list<LorentzVector> Jets;
                const Jets& jets = jetpro.getLorentzJets();
	        #endif    

		///remove leptons from list of jets
		LorentzVector jetaxis;
		_jetaxes.clear();
		for (Jets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
		  
		  if (dilep[0]->getMomentum().deltaR(*jt) > R_isol
		      && dilep[1]->getMomentum().deltaR(*jt) > R_isol) { 
		    jetaxis.setPx(jt->px());
		    jetaxis.setPy(jt->py());
		    jetaxis.setPz(jt->pz());
		    jetaxis.setE(jt->e());  
		    _jetaxes.push_back(jetaxis);
		  }
		}
		
		
		///Fill jets Et and eta distributions
		for (unsigned int i=0; i<_jetaxes.size(); ++i) {
		  if (_jetaxes[i].et()>20. && 
		      fabs(_jetaxes[i].pseudoRapidity())<1.5) {
		    _histJetsPt->fill(_jetaxes[i].et(), event.weight());
		    _histJetsEta->fill(_jetaxes[i].pseudoRapidity(), event.weight());
		  }
		}
		
		/// determine secondary vertices
		const SVertex& svtx =  event.applyProjection(_svtxproj);
		
		const vector<LorentzVector>& taggedJets = svtx.getTaggedJets();

		log << Log::DEBUG << "taggedJets.size()=" << taggedJets.size() << endl;

		///Fill tagged jets Et and eta distributions
		for (unsigned int i=0; i<taggedJets.size(); ++i) {
		  if (taggedJets[i].et()>20. && 
		      fabs(taggedJets[i].pseudoRapidity())<1.5) {
		    _histJetsPt->fill(taggedJets[i].et(), event.weight());
		    _histJetsEta->fill(taggedJets[i].pseudoRapidity(), event.weight());
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
  void HepEx0605099::finalize() { 

    normalize(_histJetsPt);
    normalize(_histJetsEta);
    normalize(_histbJetsPt);
    normalize(_histbJetsEta);

  }

}
