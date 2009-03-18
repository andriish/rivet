// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/ExampleTree.hh"


namespace Rivet {

  #ifndef HAVE_ROOT
  

  ExampleTree::ExampleTree() { }
  void ExampleTree::init() {
    getLog() << Log::WARN << "Rivet was not compiled against ROOT. ExampleTree will do nothing." << endl;
  }
  void ExampleTree::analyze(const Event& event) { }
  void ExampleTree::finalize() { }


  #else


  ExampleTree::ExampleTree() {
    const FinalState fs(-4.0, 4.0, 0.0*GeV);
    addProjection(fs, "FS");
    addProjection(ChargedLeptons(fs), "ChLeptons");
    addProjection(FastJets(fs, FastJets::KT, 0.7), "Jets");
    
    /// Veto neutrinos, antineutrinos and LSP
    VetoedFinalState vfs(fs);
    vfs
      .addVetoDetail(NU_E, 10.0*GeV, 50.0*GeV)
      .addVetoPairId(NU_MU)
      .addVetoPairId(NU_TAU)
      .addVetoId(1000022); // LSP
    addProjection(vfs, "VFS");
    addProjection(TotalVisibleMomentum(vfs), "TotalVisMom");

    ZFinder zs(fs, ELECTRON, 80*GeV, 100*GeV, 0.2);
    addProjection(zs, "Zs");
  }


  void ExampleTree::init() {
    // Choose cuts
    _jet_pt_cut = 20*GeV;
    _subj_pt_cut = 20*GeV;
    _lepton_pt_cut = 20*GeV;
    _store_partons = true;

    _treeFileName = "rivetTree.root";

    // Create a file for the Tree
    _treeFile = new TFile(_treeFileName, "recreate");

    // Book the ntuple.
    _rivetTree = new TTree("Rivet Tree", "Rivet Example Tree");

    // Event number 
    _rivetTree->Branch("nevt", &_nevt, "nevt/I");

    // Vector bosons
    _rivetTree->Branch("nvb", &_nvb, "nvb/I");
    _rivetTree->Branch("vbtype", &_vbtype, "vbtype[nvb]/I");
    _rivetTree->Branch("vbvec", &_vbvec, "vbvec[nvb][4]/F");

    _rivetTree->Branch("njet", &_njet, "njet/I");
    _rivetTree->Branch("vjet", &_vjet, "vjet[njet][4]/F");

    _rivetTree->Branch("nsub", &_nsub, "nsub/I");
    _rivetTree->Branch("sjet3", &_sjet3, "sjet3[nsub][4]/F");
    _rivetTree->Branch("ysubsj", &_ysubsj, "ysubsj[nsub][4]/F");

    _rivetTree->Branch("nlep", &_nlep, "nlep/I");
    _rivetTree->Branch("vlep", &_vlep, "vlep[nlep][4]/F");
    _rivetTree->Branch("leptype", &_leptype, "leptype[nlep][3]/F");

    _rivetTree->Branch("npart", &_npart, "npart/I");
    _rivetTree->Branch("ppart", &_ppart, "ppart[npart][4]/F");
    _rivetTree->Branch("pid", &_pid, "pid[npart]/I");
    _rivetTree->Branch("mo", &_mo, "mo[npart]/I");  // first mother.

    _rivetTree->Branch("esumr", &_esumr, "esumr[4]/F");
  }



  // Do the analysis
  void ExampleTree::analyze(const Event & event) {
    const GenEvent& ev = event.genEvent();
    _nevt = ev.event_number();

    // Get the vector bosons
    _nvb = 0;
    const FinalState& zs = applyProjection<FinalState>(event, "Zs");
    foreach (const Particle& p, zs.particles()) {
      const FourMomentum p4 = p.momentum();
      _vbvec[_nvb][0] = p4.E()/GeV;
      _vbvec[_nvb][1] = p4.px()/GeV;
      _vbvec[_nvb][2] = p4.py()/GeV;
      _vbvec[_nvb][3] = p4.pz()/GeV;
      _vbtype[_nvb]   = 1;
      ++_nvb;
    }

    // Get the partons. This is generator-dependent and should not be
    // used in normal analyses.
    _npart = 0;
    if (_store_partons) {
      for (GenEvent::particle_const_iterator pi = event.genEvent().particles_begin(); 
           pi != event.genEvent().particles_end(); ++pi ) {
        // Only include particles which are documentation line (status >1) 
        // The result/meaning will be generator dependent.
        if ( (*pi)->status() >= 2 ) {
          const FourMomentum p4 = (*pi)->momentum();
          _ppart[_npart][1] = p4.px();
          _ppart[_npart][2] = p4.py();
          _ppart[_npart][3] = p4.pz();
          _ppart[_npart][0] = p4.E();
          _pid[_npart] = (*pi)->pdg_id();
          const GenVertex* vertex = (*pi)->production_vertex();
          // Get the first mother
          if (vertex) {
            if (vertex->particles_in_size()>0) {
              GenVertex::particles_in_const_iterator p1 = vertex->particles_in_const_begin();
              _mo[_npart] = (*p1)->pdg_id();
            } else {
              _mo[_npart] = 0;
            }
          } else {
            _mo[_npart] = 0;
          }
          getLog() << Log::DEBUG << _npart << ":" << _pid[_npart] << endl;
          ++_npart;
        }
      }
    }
    
    
    // Get the jets in decreasing pT order.
    const FastJets& jets = applyProjection<FastJets>(event, "Jets");
    PseudoJets jetList = jets.pseudoJetsByPt();
    _njet = 0;
    _nsub = 0;
    foreach (const fastjet::PseudoJet& j, jetList) {
      if (j.perp() > _jet_pt_cut) {
        _vjet[_njet][0] = j.e()/GeV;
        _vjet[_njet][1] = j.px()/GeV;
        _vjet[_njet][2] = j.py()/GeV;
        _vjet[_njet][3] = j.pz()/GeV;
        if (j.perp() > _subj_pt_cut) {
          _sjet3[_nsub][0] = j.e()/GeV;
          _sjet3[_nsub][1] = j.px()/GeV;
          _sjet3[_nsub][2] = j.py()/GeV;
          _sjet3[_nsub][3] = j.pz()/GeV;
          const vector<double> ys = jets.ySubJet(j);
          for (size_t i = 0; i < 4; ++i){
            if (ys.size() > i) {
              _ysubsj[_nsub][i] = ys.at(i);
            } else {
              _ysubsj[_nsub][i] = 0;
            }
          }
          ++_nsub;	 
        }
        ++_njet;
      }
    }
    
    // Loop over leptons
    _nlep = 0;
    const ChargedLeptons& cl = applyProjection<ChargedLeptons>(event, "ChLeptons");
    foreach (const Particle& p, cl.chargedLeptons()) {
      const FourMomentum p4 = p.momentum();
      if (p4.pT() > _lepton_pt_cut) {
		_vlep[_nlep][0] = p4.E()/GeV;
        _vlep[_nlep][1] = p4.px()/GeV;
        _vlep[_nlep][2] = p4.py()/GeV;
		_vlep[_nlep][3] = p4.pz()/GeV;
        ++_nlep;
      }
    }
    
    // Missing Et/total energy
    const TotalVisibleMomentum& tvm = applyProjection<TotalVisibleMomentum>(event, "TotalVisMom");
	_esumr[0] = tvm.momentum().E()/GeV;
    _esumr[1] = tvm.momentum().px()/GeV;
    _esumr[2] = tvm.momentum().py()/GeV;
	_esumr[3] = tvm.momentum().pz()/GeV;

    // Finally fill the tree
    _rivetTree->Fill();
  }
  
  
  // Finalize
  void ExampleTree::finalize() { 
    // Write the tree to file.
    _rivetTree->Write();
  }

  
  #endif
  
}
