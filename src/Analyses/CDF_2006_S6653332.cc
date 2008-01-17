// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2006_S6653332.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"


namespace Rivet {


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
    if (fabs(pv.getPVPosition().z()) < _pvzmax) {

    const TotalVisibleMomentum& caloMissEt = event.applyProjection(_calmetproj);
    log << Log::DEBUG << "CaloMissEt.getMomentum().pT() = " << caloMissEt.getMomentum().pT() << endl;

    if (caloMissEt.getMomentum().pT() < _metmax && caloMissEt.getSET() < _setmax) {
      // Get leptons in event
      const ChargedLeptons& chlep = event.applyProjection(_chleproj);
      const ParticleVector& leps = chlep.chargedLeptons();
      
      // Sort 2 leading leptons in pT
      vector<ParticleVector::const_iterator> dilep;
      for (ParticleVector::const_iterator p = leps.begin(); p != leps.end(); ++p) {
        if (dilep.empty()) {
          dilep.push_back(p);
        } 
        else if (dilep.size() == 1) {
          if (p->getMomentum().pT() > dilep[0]->getMomentum().pT()) { 
            dilep.push_back(dilep[0]);
            dilep[0]= p;
          } else {
            dilep.push_back(p);
          }
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
      if (dilep.size() == 2 && dilep[1]->getMomentum().pT() > _lepPtmin 
          && dilep[0]->getPdgId() == -dilep[1]->getPdgId() ) {

        // Muon and electron geometrical acceptance requirement
        const double eta0 = dilep[0]->getMomentum().pseudorapidity();
        const double eta1 = dilep[1]->getMomentum().pseudorapidity();
        if ((abs(dilep[0]->getPdgId() == ELECTRON) && fabs(eta0) < _eletamax && fabs(eta1) < _eletamax) ||
            (abs(dilep[0]->getPdgId() == MUON) && fabs(eta0) < _muetamax && fabs(eta1) < _muetamax)) {

          // Dilepton invariant mass cut
          FourMomentum ll = dilep[0]->getMomentum() + dilep[1]->getMomentum();
          log << Log::DEBUG << "ll.mass()=" << ll.mass() << endl;
          if (_mllmin < ll.mass()  &&  ll.mass() < _mllmax) {
	    
            log << Log::DEBUG << "Invariant di-lepton mass ll.mass() = " << ll.mass() << endl;
            
            double rapmax;
            if (abs(dilep[0]->getPdgId()) == ELECTRON) rapmax = _trigeletamax; //electron
            else rapmax = _trigmuetamax; //muon
            
            // More restrictive pT requirement on (one) trigger lepton
            if ( (dilep[0]->getMomentum().pT() > _triglepPtmin && 
                  fabs(dilep[0]->getMomentum().pseudorapidity()) < rapmax) || 
                 (dilep[1]->getMomentum().pT() > _triglepPtmin && 
                  fabs(dilep[1]->getMomentum().pseudorapidity()) < rapmax)) {
              
              
              const VetoedFinalState& vfs = event.applyProjection(_vfsproj);
              // Check if leptons are isolated
              vector<double> Ehad(2, 0.); //dim=2 for 2 leptons
              vector<double> Eem(2, 0.); //dim=2 for 2 leptons

              for (ParticleVector::const_iterator p = vfs.particles().begin(); p != vfs.particles().end(); ++p) {
                const int pid = p->getPdgId();
                const int abspid = abs(pid);
                // Determine energy depositions around lepton(s)
                for (size_t lind = 0; lind < 2; ++lind) {
                  if (deltaR(p->getMomentum(), dilep[lind]->getMomentum()) < _Rlepisol) {
                    if (abspid == ELECTRON || abspid == MUON || abspid == PHOTON) {
                      Eem[lind] += p->getMomentum().E();
                    } else if (PID::threeCharge(pid) != 0) { 
                      // Charged hadron: had + em calo fraction
                      Eem[lind] += (1.-_fh)*p->getMomentum().E();
                      Ehad[lind] += _fh*p->getMomentum().E();
                    } else { 
                      // Neutral hadron assumed to have vanishing EM fraction
                      Ehad[lind] += p->getMomentum().E();
                    }
                  }
                }
              }
              
              bool l_isol = true; // Both leptons have to be isolated 
              // => one veto variable suffices
              for (size_t lind = 0; lind < 2; ++lind) {
                if (abs(dilep[lind]->getPdgId()) == ELECTRON) { // electron
                  if (!(Eem[lind]>0.)) l_isol = false;
                  else if (Ehad[lind]/Eem[lind] > _fhfemconst + _fhfemslope*dilep[lind]->getMomentum().E()) l_isol = false;
                }
                else { // muon
                  if (dilep[lind]->getMomentum().E() <= _muEsep) {
                    if (Ehad[lind]-dilep[lind]->getMomentum().E() > _muEhMin
                        || Eem[lind]-dilep[lind]->getMomentum().E() > _muEemMin) l_isol = false;
                  }
                  else { // p > 100.
                    if (Ehad[lind]-dilep[lind]->getMomentum().E() > // E_had = 0.0280 
                        _muEhMin+(dilep[lind]->getMomentum().E()-_muEsep)*_muEhslope
                        || Eem[lind]-dilep[lind]->getMomentum().E() > // E_em = 0.0115 
                        _muEemMin+(dilep[lind]->getMomentum().E()-_muEsep)*_muEemslope) l_isol = false;
                  }
                }
              }

              if (l_isol) { // Both leptons are isolated
                log << Log::DEBUG << "Both leptons are isolated" << endl;
                const FastJets& jetpro = event.applyProjection(_jetsproj);
                const PseudoJets& jets = jetpro.getPseudoJetsPt();
                
                // Remove leptons from list of jets
                _jetaxes.clear();
                for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
                  const FourMomentum j(jt->E(), jt->px(), jt->py(), jt->pz());
                  if (deltaR(dilep[0]->getMomentum(), j) > _Rlepisol && 
                      deltaR(dilep[1]->getMomentum(), j) > _Rlepisol) { _jetaxes.push_back(j); }
                }
                
                // Fill jets Et and eta distributions.
                for (size_t i=0; i < _jetaxes.size(); ++i) {
                  if (_jetaxes[i].Et() > _jetPtmin && 
                      fabs(_jetaxes[i].pseudorapidity()) < _jetetamax) {
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
                  if (taggedJets[i].Et() > _jetPtmin && 
                      fabs(taggedJets[i].pseudorapidity()) < _jetetamax) {
                    _histJetsPt->fill(taggedJets[i].Et(), event.weight());
                    _histJetsEta->fill(taggedJets[i].pseudorapidity(), event.weight());
                  }
                }
                
              } // Both leptons are isolated
            } // Dilepton invariant mass cut
          } // mu,e geometrical acceptance requirement
        } // At least one trigger lepton (eta, pT) requirement	
      } // Opposite charged leptons, pT>10
      
    } // Cal. missing Et and scalar Et (Ht) cut
    
    } // Primary vertex z-position requirement
    
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
