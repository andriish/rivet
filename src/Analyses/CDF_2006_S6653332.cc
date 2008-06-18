// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2006_S6653332.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"


namespace Rivet {


  // Book histograms
  void CDF_2006_S6653332::init() {
    _histJetsPt   = bookHistogram1D(4, 1, 1, "All jets Pt");
    _histJetsEta  = bookHistogram1D(5, 1, 1, "All jets Pseudorapidity");
    _histbJetsPt  = bookHistogram1D(6, 1, 1, "b jets Pt");
    _histbJetsEta = bookHistogram1D(7, 1, 1, "b jets Pseudorapidity");
  }


  // Do the analysis
  void CDF_2006_S6653332::analyze(const Event& event) {
    Log log = getLog();
    
    // Find primary vertex 
    const PVertex& pv = applyProjection<PVertex>(event, "PV");
    if (fabs(pv.getPVPosition().z())/mm > _pvzmax) {
      vetoEvent(event);
    }

    // Veto on missing Et    
    const TotalVisibleMomentum& caloMissEt = applyProjection<TotalVisibleMomentum>(event, "CalMET");
    log << Log::DEBUG << "Missing pT = " << caloMissEt.getMomentum().pT() << endl;
    if (caloMissEt.getMomentum().pT() > _metmax || caloMissEt.getSET() > _setmax) {
      vetoEvent(event);
    }

    // Get leptons in event
    const ParticleVector& leps = applyProjection<ChargedLeptons>(event, "ChLeptons").chargedLeptons();
    if (leps.empty()) vetoEvent(event);

    // Sort 2 leading leptons in pT 
    vector<ParticleVector::const_iterator> dilep;
    for (ParticleVector::const_iterator p = leps.begin(); p != leps.end(); ++p) {
      if (dilep.empty()) {
        dilep.push_back(p);
      } else if (dilep.size() == 1) {
        if (p->getMomentum().pT() > dilep[0]->getMomentum().pT()) { 
          dilep.push_back(dilep[0]);
          dilep[0] = p;
        } else {
          dilep.push_back(p);
        }
      }
      else { // dilep.size == 2
        if (p->getMomentum().pT() > dilep[0]->getMomentum().pT() )  { 
          dilep[1] = dilep[0];
          dilep[0] = p;
        } else if (p->getMomentum().pT() > dilep[1]->getMomentum().pT() )  { 
          dilep[1] = p; 
        }
      }
    }
    
    // Require exactly 2 leptons, same flavour, opposite charge, pTmin > 10 GeV.
    if (dilep.size() != 2 || 
        dilep[1]->getMomentum().pT() < _lepPtmin ||
        dilep[0]->getPdgId() != -dilep[1]->getPdgId() ) {
      vetoEvent(event);
    }

    // Muon and electron geometrical acceptance requirement
    const long dilepId = dilep[0]->getPdgId();
    const double eta0 = dilep[0]->getMomentum().pseudorapidity();
    const double eta1 = dilep[1]->getMomentum().pseudorapidity();
    if (abs(dilepId) == ELECTRON) {
      if (fabs(eta0) > _eletamax || fabs(eta1) > _eletamax) vetoEvent(event);
    } else if (abs(dilepId) == MUON) {
      if (fabs(eta0) > _muetamax || fabs(eta1) > _muetamax) vetoEvent(event);
    } else {
      vetoEvent(event);
    }
    // This was the buggy former version, which takes abs(bool) and hence 
    // probably misses positrons and antimuons:
    // if (( abs(dilep[0]->getPdgId() == ELECTRON) && fabs(eta0) < _eletamax && fabs(eta1) < _eletamax ) ||
    //     (abs(dilep[0]->getPdgId() == MUON) && fabs(eta0) < _muetamax && fabs(eta1) < _muetamax)) {
        
    // Dilepton invariant mass cut
    FourMomentum ll = dilep[0]->getMomentum() + dilep[1]->getMomentum();
    log << Log::DEBUG << "Dilepton mass = " << ll.mass() << endl;
    if (ll.mass() < _mllmin || ll.mass() > _mllmax) {
      vetoEvent(event);
    }
    
    // More restrictive pT requirement on (one) trigger lepton
    double rapmax = _trigeletamax; //  electron value by default
    if (abs(dilepId) == MUON) rapmax = _trigmuetamax;
    const FourMomentum& mom0 = dilep[0]->getMomentum();
    const FourMomentum& mom1 = dilep[1]->getMomentum();
    if ( (mom0.pT() < _triglepPtmin || fabs(mom0.pseudorapidity()) > rapmax) && 
         (mom1.pT() < _triglepPtmin || fabs(mom1.pseudorapidity()) > rapmax)) {
      vetoEvent(event);
    }

    // Check if leptons are isolated
    const VetoedFinalState& vfs = applyProjection<VetoedFinalState>(event, "VFS");
    vector<double> Ehad(2, 0.0), Eem(2, 0.0); // dim=2 for 2 leptons
    for (ParticleVector::const_iterator p = vfs.particles().begin(); p != vfs.particles().end(); ++p) {
      const long pid = p->getPdgId();
      const unsigned long abspid = abs(pid);
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
      } else { // muon
        if (dilep[lind]->getMomentum().E() <= _muEsep) {
          if (Ehad[lind]-dilep[lind]->getMomentum().E() > _muEhMin
              || Eem[lind]-dilep[lind]->getMomentum().E() > _muEemMin) l_isol = false;
        } else { // p > 100.
          if (Ehad[lind]-dilep[lind]->getMomentum().E() > // E_had = 0.0280 
              _muEhMin+(dilep[lind]->getMomentum().E()-_muEsep)*_muEhslope
              || Eem[lind]-dilep[lind]->getMomentum().E() > // E_em = 0.0115 
              _muEemMin+(dilep[lind]->getMomentum().E()-_muEsep)*_muEemslope) l_isol = false;
        }
      }
    }
    if (!l_isol) vetoEvent(event);
    log << Log::DEBUG << "Both leptons are isolated" << endl;

    // Get list of jets and remove leptons
    const PseudoJets& jets = applyProjection<FastJets>(event, "Jets").getPseudoJetsByPt();
    _jetaxes.clear();
    for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      const FourMomentum j(jt->E(), jt->px(), jt->py(), jt->pz());
      if (deltaR(dilep[0]->getMomentum(), j) > _Rlepisol && 
          deltaR(dilep[1]->getMomentum(), j) > _Rlepisol) { 
        _jetaxes.push_back(j); 
      }
    }
      
    // Fill jets Et and eta distributions.
    for (size_t i = 0; i < _jetaxes.size(); ++i) {
      if (_jetaxes[i].Et() > _jetPtmin && 
          fabs(_jetaxes[i].pseudorapidity()) < _jetetamax) {
        _histJetsPt->fill(_jetaxes[i].Et(), event.weight());
        _histJetsEta->fill(_jetaxes[i].pseudorapidity(), event.weight());
      }
    }
      
    // Determine secondary vertices
    const SVertex& svtx = applyProjection<SVertex>(event, "SVtx");
    const vector<FourMomentum>& taggedJets = svtx.getTaggedJets();
    log << Log::DEBUG << "taggedJets.size()=" << taggedJets.size() << endl;
      
    // Fill tagged jets Et and eta distributions
    for (vector<FourMomentum>::const_iterator j = taggedJets.begin(); j != taggedJets.end(); ++j) {
      if (j->Et() > _jetPtmin && fabs(j->pseudorapidity()) < _jetetamax) {
        _histJetsPt->fill(j->Et(), event.weight());
          _histJetsEta->fill(j->pseudorapidity(), event.weight());
      }
    }
  }
    
    
  
  // Finalize
  void CDF_2006_S6653332::finalize() { 
    normalize(_histJetsPt);
    normalize(_histJetsEta);
    normalize(_histbJetsPt);
    normalize(_histbJetsEta);
  }
  
  
}
