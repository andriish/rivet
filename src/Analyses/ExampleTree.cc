// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/ExampleTree.hh"

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


namespace Rivet {

  #ifndef HAVE_ROOT
  void ExampleTree::init() {
    getLog() << Log::WARN << "Rivet was not compiled against ROOT. ExampleTree will do nothing." << endl;
  }
  void ExampleTree::analyze(const Event & event) { }
  void ExampleTree::finalize() { }
  #else


  void ExampleTree::init() {
    // Book histograms
    _jet_pt_cut = 20;
    _subj_pt_cut = 20;
    _lepton_pt_cut = 20;
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
    Log log = getLog();
    log << Log::DEBUG << "Filling the ntuple" << endl;

    const GenEvent& ev = event.genEvent();

    // Event number
    _nevt = ev.event_number();

    // Jets.
    const FastJets& jets = event.applyProjection(_jetsproj);

    // Leptons
    const ChargedLeptons& cl = event.applyProjection(_chgleptonsproj);

    // Missing Et/total energy
    const TotalVisibleMomentum& tvm = event.applyProjection(_totvismomproj);

    // Vector bosons.
    const WZandh& wzh = event.applyProjection(_wzandhproj);

    // Get the vector bosons
    _nvb = 0;
    for (ParticleVector::const_iterator p = wzh.Zees().begin(); p != wzh.Zees().end(); ++p) {
      const FourMomentum p4 = p->getMomentum();
      _vbvec[_nvb][1] = p4.px();
      _vbvec[_nvb][2] = p4.py();
      _vbvec[_nvb][3] = p4.pz();
      _vbvec[_nvb][0] = p4.E();
      _vbtype[_nvb]   = 1;
      ++_nvb;
    }

    // Get the partons. This is generator dependent and should not be
    // used in normal analyses.
    _npart = 0;
    if (_store_partons) {

      for (GenEvent::particle_const_iterator pi = event.genEvent().particles_begin(); 
           pi != event.genEvent().particles_end(); ++pi ) {
        // Only include particles which are documentation line (status = 3) 
        // The result/meaning will be generator dependent.
        if ( (*pi)->status() >= 2 ) {
          const FourMomentum p4 = (*pi)->momentum();
          _ppart[_npart][1] = p4.px();
          _ppart[_npart][2] = p4.py();
          _ppart[_npart][3] = p4.pz();
          _ppart[_npart][0] = p4.E();
          HepPDT::ParticleID pInfo = (*pi)->pdg_id();
          _pid[_npart] = pInfo.pid();
          const GenVertex* vertex = (*pi)->production_vertex();
          // get the first mother
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
          
          //const GenParticle mother = **(vertex->particles_in_const_begin());
          log << Log::DEBUG << _npart << ":" << _pid[_npart] << endl;
          ++_npart;
        }
      }
    }
    
    
    // Get the jets in decreasing ET order.
    typedef vector<fastjet::PseudoJet> Jets;
    Jets jetList = jets.getJets();
    _njet = 0;
    _nsub = 0;
    for (Jets::const_iterator j = jetList.begin(); j != jetList.end(); ++j) {
      if (j->perp() > _jet_pt_cut) {
        _vjet[_njet][1] = j->px();
        _vjet[_njet][2] = j->py();
        _vjet[_njet][3] = j->pz();
        _vjet[_njet][0] = j->e();
        if (j->perp() > _subj_pt_cut) {
          _sjet3[_nsub][1] = j->px();
          _sjet3[_nsub][2] = j->py();
          _sjet3[_nsub][3] = j->pz();
          _sjet3[_nsub][0] = j->e();
          vector<double> ys = jets.getYSubJet(*j);
          for (size_t i=0; i<5; ++i){
            if (ys.size()>i) {
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
    for (ParticleVector::const_iterator p = cl.chargedLeptons().begin(); p != cl.chargedLeptons().end(); ++p) {
      const FourMomentum p4 = p->getMomentum();
      if (p4.pT() > _lepton_pt_cut) {
        _vlep[_nlep][1] = p4.px();
        _vlep[_nlep][2] = p4.py();
		_vlep[_nlep][3] = p4.pz();
		_vlep[_nlep][0] = p4.E();
        ++_nlep;
      }
    }
    
    // Total/missing energy.  
    _esumr[1] = tvm.getMomentum().px();
    _esumr[2] = tvm.getMomentum().py();
    _esumr[3] = tvm.getMomentum().pz();
    _esumr[0] = tvm.getMomentum().E();
    
    // Finished...
    log << Log::DEBUG << "Finished analyzing" << endl;
    
    _rivetTree->Fill();
  }
  
  
  // Finalize
  void ExampleTree::finalize() { 
    // Write the tree to file.
    _rivetTree->Write();
  }
  
  #endif
  
}
