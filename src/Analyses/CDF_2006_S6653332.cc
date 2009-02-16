// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2006_S6653332.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"


namespace Rivet {

  // Constructor
  // Analysis cut values:
  //  - _pvzmax: cut on primary vertex \f$ z \f$ position ( \f$ z(\text{PV}) < 60 \text{cm} \f$ );
  //  - _met: cut on missing transverse energy ( \f$ MET < 25 \text{GeV} \f$ )
  //  - _set: cut on scalar sum of transverse energy ( \f$ SET < 150 \text{GeV} \f$ )
  //  - _lepPtmin: cut on leptons minimum transverse momentum ( \f$ pT(\text{lep.}) > 10 \text{GeV} \f$ )
  //  - _eletamax: cut on max. pseudorapidity of electrons/calorimeter ( \f$ |eta(el.)| < 3.6 \f$ ) 
  //  - _muetamax: cut on max. pseudorapidity of muons/muon chambers ( \f$ |eta(el.)| < 1.5 \f$ ) 
  //  - _mllmin: cut on min invariant dilepton mass ( \f$ m(ll) > 66 \text{GeV} \f$ )
  //  - _mllmax: cut on max invariant dilepton mass ( \f$ m(ll) <116 \text{GeV} \f$ )
  //  - _triglepPtmin: cut on trigger lepton min. transverse momentum ( \f$ pT > 18 \text{GeV} \f$ )
  //  - _triglepetamax: cut on trigger lepton max pseudorapidity ( \f$ |eta| < 1.1 \f$ )
  //  - _Rlepisol: cut on lepton isolation radius ( \f$ R < 0.4 \f$ )
  //  - hadronic calorimeter fraction of charged hadrons, needed for cuts on hadr. + el. mag. fractions
  //  - hadronic/el.magn. calo. fraction cut const ( \f$ fhfemconst = 0.055 \f$ )
  //  - hadronic/el.magn. calo. fraction cut slope ( \f$ fhfemslope = 0.00045 \f$ )
  //  - muon Energy cut case separation ( \f$ muEsep = 100 \text{GeV} \f$ )
  //  - muon min. E hadr.isol. threshold ( \f$ muEhMin = 6 \text{GeV} \f$ )
  //  - muon min. E el.mag. isol. threshold ( \f$ muEemMin = 2 \text{GeV} \f$ )
  //  - muon had. Energy isol. cut slope ( \f$ muEhslope = 0.0280 \f$ )
  //  - muon el.mag. Energy isol. cut slope ( \f$ muEemslope = 0.0115 \f$ )
  //  - jets transverse momentum cut ( \f$ p_\text{T,jet} > 20 \text{GeV} \f$ )
  //  - jets pseudorapidity cut ( \f$ \eta_\text{jet} < 1.5 \f$ )
  //  - _vtxjetRmax: Max. distance in \f$ (\eta, \phi) \f$ space between vertex vis. momentum and jet to be probed ( \f$ R_\text{vtx-jet} < 0.7 \f$ )  
  //  - _trketamax: Tracker geometrical acceptance ( \f$ \eta < 2.0 \f$ )
  //  - _ipres: Impact Parameter resolution ( \f$ \Delta{\text{IP}} = 34e-3 \text{mm} \f$ ), including beam spot
  //  - _dlsmin: cut on Decay Length Significance ( \f$ l/\Delta{l} > 7.5 \f$ )
  //  - _dlsres: Decay Length Significance resolution (assumed to be ( \f$ \Delta{l} = 34e-3 \text{mm} \f$ ))
  CDF_2006_S6653332::CDF_2006_S6653332()
    : _pvzmax(600*mm), _metmax(25*GeV), _setmax(150*GeV), _lepPtmin(10*GeV),
      _eletamax(3.5), _muetamax(1.5), _mllmin(66*GeV), _mllmax(116*GeV),
      _triglepPtmin(18*GeV), _trigeletamax(1.1), _trigmuetamax(1.0),
      _Rlepisol(0.4), _fh(0.8), _fhfemconst(0.055), _fhfemslope(0.00045),
      _muEsep(100*GeV), _muEhMin(6*GeV), _muEemMin(2*GeV), _muEhslope(0.0280),
      _muEemslope(0.0115), _jetPtmin(20*GeV), _jetetamax(1.5),
      _vtxjetRmax(0.7), _trketamax(2.0), _ipres(34e-3*mm), _dlsmin(7.5), _dlsres(34e-3*mm)
  {
    setBeams(PROTON, ANTIPROTON);

    // Veto (anti)neutrinos, and muons with \f$ p_T \f$ above 1.0 GeV
    /// @todo If we allow VFS constructor to specify eta and pT ranges, we can
    /// bypass making this FS
    FinalState fs(-3.6, 3.6);
    addProjection(fs, "FS");
    VetoedFinalState vfs(fs);
    vfs
      .addVetoPairId(NU_E)
      .addVetoPairId(NU_MU)
      .addVetoPairId(NU_TAU)
      .addVetoDetail(MUON, 1.0*GeV, MAXDOUBLE);
    addProjection(vfs, "VFS");
    addProjection(FastJets(vfs), "Jets");
    addProjection(ChargedFinalState(vfs), "ChFS");
    addProjection(TotalVisibleMomentum(vfs), "CalMET");
    addProjection(ChargedLeptons(vfs), "ChLeptons");
    addProjection(PVertex(), "PV");
    addProjection(SVertex(getProjection<ChargedFinalState>("ChFS"),
                          _jetaxes, _vtxjetRmax, _trketamax, 
                          _ipres, _dlsmin, _dlsres), "SVtx");
  }



  // Book histograms
  void CDF_2006_S6653332::init() {
    _histJetsPt   = bookHistogram1D(4, 1, 1, "All jets $p_\\perp$");
    _histJetsEta  = bookHistogram1D(5, 1, 1, "All jets pseudorapidity, $\\eta$");
    _histbJetsPt  = bookHistogram1D(6, 1, 1, "$b$ jets $p_\\perp$");
    _histbJetsEta = bookHistogram1D(7, 1, 1, "$b$ jets pseudorapidity, $\\eta$");
  }


  // Do the analysis
  void CDF_2006_S6653332::analyze(const Event& event) {
    Log log = getLog();
    
    // Find primary vertex 
    const PVertex& pv = applyProjection<PVertex>(event, "PV");
    if (fabs(pv.position().z())/mm > _pvzmax) {
      vetoEvent(event);
    }

    // Veto on missing Et    
    const TotalVisibleMomentum& caloMissEt = applyProjection<TotalVisibleMomentum>(event, "CalMET");
    log << Log::DEBUG << "Missing pT = " << caloMissEt.momentum().pT() << endl;
    if (caloMissEt.momentum().pT() > _metmax || caloMissEt.scalarET() > _setmax) {
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
        if (p->momentum().pT() > dilep[0]->momentum().pT()) { 
          dilep.push_back(dilep[0]);
          dilep[0] = p;
        } else {
          dilep.push_back(p);
        }
      }
      else { // dilep.size == 2
        if (p->momentum().pT() > dilep[0]->momentum().pT() )  { 
          dilep[1] = dilep[0];
          dilep[0] = p;
        } else if (p->momentum().pT() > dilep[1]->momentum().pT() )  { 
          dilep[1] = p; 
        }
      }
    }
    
    // Require exactly 2 leptons, same flavour, opposite charge, pTmin > 10 GeV.
    if (dilep.size() != 2 || 
        dilep[1]->momentum().pT() < _lepPtmin ||
        dilep[0]->pdgId() != -dilep[1]->pdgId() ) {
      vetoEvent(event);
    }

    // Muon and electron geometrical acceptance requirement
    const long dilepId = dilep[0]->pdgId();
    const double eta0 = dilep[0]->momentum().pseudorapidity();
    const double eta1 = dilep[1]->momentum().pseudorapidity();
    if (abs(dilepId) == ELECTRON) {
      if (fabs(eta0) > _eletamax || fabs(eta1) > _eletamax) vetoEvent(event);
    } else if (abs(dilepId) == MUON) {
      if (fabs(eta0) > _muetamax || fabs(eta1) > _muetamax) vetoEvent(event);
    } else {
      vetoEvent(event);
    }
    // This was the buggy former version, which takes abs(bool) and hence 
    // probably misses positrons and antimuons:
    // if (( abs(dilep[0]->pdgId() == ELECTRON) && fabs(eta0) < _eletamax && fabs(eta1) < _eletamax ) ||
    //     (abs(dilep[0]->pdgId() == MUON) && fabs(eta0) < _muetamax && fabs(eta1) < _muetamax)) {
        
    // Dilepton invariant mass cut
    FourMomentum ll = dilep[0]->momentum() + dilep[1]->momentum();
    log << Log::DEBUG << "Dilepton mass = " << ll.mass() << endl;
    if (ll.mass() < _mllmin || ll.mass() > _mllmax) {
      vetoEvent(event);
    }
    
    // More restrictive pT requirement on (one) trigger lepton
    double rapmax = _trigeletamax; //  electron value by default
    if (abs(dilepId) == MUON) rapmax = _trigmuetamax;
    const FourMomentum& mom0 = dilep[0]->momentum();
    const FourMomentum& mom1 = dilep[1]->momentum();
    if ( (mom0.pT() < _triglepPtmin || fabs(mom0.pseudorapidity()) > rapmax) && 
         (mom1.pT() < _triglepPtmin || fabs(mom1.pseudorapidity()) > rapmax)) {
      vetoEvent(event);
    }

    // Check if leptons are isolated
    const VetoedFinalState& vfs = applyProjection<VetoedFinalState>(event, "VFS");
    vector<double> Ehad(2, 0.0), Eem(2, 0.0); // dim=2 for 2 leptons
    for (ParticleVector::const_iterator p = vfs.particles().begin(); p != vfs.particles().end(); ++p) {
      const long pid = p->pdgId();
      const unsigned long abspid = abs(pid);
      // Determine energy depositions around lepton(s)
      for (size_t lind = 0; lind < 2; ++lind) {
        if (deltaR(p->momentum(), dilep[lind]->momentum()) < _Rlepisol) {
          if (abspid == ELECTRON || abspid == MUON || abspid == PHOTON) {
            Eem[lind] += p->momentum().E();
          } else if (PID::threeCharge(pid) != 0) { 
            // Charged hadron: had + em calo fraction
            Eem[lind] += (1.-_fh)*p->momentum().E();
            Ehad[lind] += _fh*p->momentum().E();
          } else { 
            // Neutral hadron assumed to have vanishing EM fraction
            Ehad[lind] += p->momentum().E();
          }
        }
      }
    }
    
    bool l_isol = true; // Both leptons have to be isolated 
    // => one veto variable suffices
    for (size_t lind = 0; lind < 2; ++lind) {
      if (abs(dilep[lind]->pdgId()) == ELECTRON) { // electron
        if (!(Eem[lind]>0.)) l_isol = false;
        else if (Ehad[lind]/Eem[lind] > _fhfemconst + _fhfemslope*dilep[lind]->momentum().E()) l_isol = false;
      } else { // muon
        if (dilep[lind]->momentum().E() <= _muEsep) {
          if (Ehad[lind]-dilep[lind]->momentum().E() > _muEhMin
              || Eem[lind]-dilep[lind]->momentum().E() > _muEemMin) l_isol = false;
        } else { // p > 100.
          if (Ehad[lind]-dilep[lind]->momentum().E() > // E_had = 0.0280 
              _muEhMin+(dilep[lind]->momentum().E()-_muEsep)*_muEhslope
              || Eem[lind]-dilep[lind]->momentum().E() > // E_em = 0.0115 
              _muEemMin+(dilep[lind]->momentum().E()-_muEsep)*_muEemslope) l_isol = false;
        }
      }
    }
    if (!l_isol) vetoEvent(event);
    log << Log::DEBUG << "Both leptons are isolated" << endl;

    // Get list of jets and remove leptons
    const PseudoJets& jets = applyProjection<FastJets>(event, "Jets").pseudoJetsByPt();
    _jetaxes.clear();
    for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      const FourMomentum j(jt->E(), jt->px(), jt->py(), jt->pz());
      if (deltaR(dilep[0]->momentum(), j) > _Rlepisol && 
          deltaR(dilep[1]->momentum(), j) > _Rlepisol) { 
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
