// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/ExampleTree.hh"
#include "Rivet/RivetAIDA.hh"
using namespace Rivet;

#include "AIDA/IHistogram1D.h"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;



/////////////////////////////////////////////////


// Book histograms
void ExampleTree::init() {

#ifndef HAVE_ROOT
  Log log = getLog();
  log << Log::WARN << "Rivet was not compiled against ROOT. ExampleTree will do nothing." << std::endl;
#else

  _jet_pt_cut = 20;
  _subj_pt_cut = 20;
  _lepton_pt_cut = 20;
  _store_partons = true;

  treeFileName = "rivetTree.root";

  // Create a file for the Tree
  treeFile = new TFile(treeFileName,"recreate");

  // book the ntuple.
  rivetTree = new TTree("Rivet Tree","Rivet Example Tree");

  // event number 
  rivetTree->Branch("nevt",&nevt,"nevt/I");

  // Vector bosons
  rivetTree->Branch("nvb",&nvb,"nvb/I");
  rivetTree->Branch("vbtype",&vbtype,"vbtype[nvb]/I");
  rivetTree->Branch("vbvec",&vbvec,"vbvec[nvb][4]/F");

  rivetTree->Branch("njet",&njet,"njet/I");
  rivetTree->Branch("vjet",&vjet,"vjet[njet][4]/F");

  rivetTree->Branch("nsub",&nsub,"nsub/I");
  rivetTree->Branch("sjet3",&sjet3,"sjet3[nsub][4]/F");
  rivetTree->Branch("ysubsj",&ysubsj,"ysubsj[nsub][4]/F");

  rivetTree->Branch("nlep",&nlep,"nlep/I");
  rivetTree->Branch("vlep",&vlep,"vlep[nlep][4]/F");
  rivetTree->Branch("leptype",&leptype,"leptype[nlep][3]/F");

  rivetTree->Branch("npart",&npart,"npart/I");
  rivetTree->Branch("ppart",&ppart,"ppart[njet][4]/F");
  rivetTree->Branch("pid",&pid,"pid[npart]/I");
  rivetTree->Branch("mo",&mo,"mo[npart]/I");  // first mother.

  rivetTree->Branch("esumr",&esumr,"esumr[4]/F");

#endif
}


// Do the analysis
void ExampleTree::analyze(const Event & event) {
#ifdef HAVE_ROOT
  Log log = getLog();
  log << Log::DEBUG << "Filling the ntuple" << endl;

  GenEvent ev = event.genEvent();

  // Event number
  nevt = ev.event_number();

  // Jets.
  const KtJets& jets = event.applyProjection(p_ktjets);

  // Leptons
  const ChargedLeptons& cl = event.applyProjection(p_chargedleptons);

  // Missing Et/total energy
  const TotalVisibleMomentum& tvm = event.applyProjection(*p_totalvisiblemomentum);

  // Vector bosons.
  const WZandh& wzh = event.applyProjection(p_wzandh);

  // Get the vector bosons
  nvb = 0;
  for (ParticleVector::const_iterator p = wzh.Zees().begin(); p != wzh.Zees().end(); ++p) {
    vbvec[nvb][1] = p->getMomentum().px();
    vbvec[nvb][2] = p->getMomentum().py();
    vbvec[nvb][3] = p->getMomentum().pz();
    vbvec[nvb][0] = p->getMomentum().e();
    vbtype[nvb]   = 1;
    nvb++;
  }

  // Get the partons. This is generator dependent and should not be
  // used in normal analyses.
  npart = 0;
  if (_store_partons) {

    for ( GenEvent::particle_const_iterator pi = event.genEvent().particles_begin(); 
	  pi != event.genEvent().particles_end(); ++pi ) {
      // Only include particles which are documentation line (status = 3) 
      // The result/meaning will be generator dependent.
      if ( (*pi)->status() >= 2 ) {
	ppart[npart][1] = (*pi)->momentum().px();
	ppart[npart][2] = (*pi)->momentum().py();
	ppart[npart][3] = (*pi)->momentum().pz();
	ppart[npart][0] = (*pi)->momentum().e();
	HepPDT::ParticleID pInfo = (*pi)->pdg_id();
	pid[npart]=pInfo.pid();
	const GenVertex *vertex = (*pi)->production_vertex();
	// get the first mother
	if (vertex) {
	  if (vertex->particles_in_size()>0) {
	    GenVertex::particles_in_const_iterator p1 = vertex->particles_in_const_begin();
	    mo[npart] = (*p1)->pdg_id();
	  } else {
	    mo[npart] = 0;
	  }
	} else {
	  mo[npart] = 0;
	}


	//const GenParticle mother = **(vertex->particles_in_const_begin());
	
	log << Log::DEBUG << npart << ":" << pid[npart] << endl;
	npart++;
      }
    }
  }


  // Get the jets in decreasing ET order.
  vector<KtJet::KtLorentzVector> jetList = jets.getJetsEt();
  njet = 0;
  nsub = 0;
  for (vector<KtJet::KtLorentzVector>::iterator j = jetList.begin(); j != jetList.end(); ++j) {
    if (j->perp()>_jet_pt_cut) {
      vjet[njet][1] = j->px();
      vjet[njet][2] = j->py();
      vjet[njet][3] = j->pz();
      vjet[njet][0] = j->e();
      if (j->perp()>_subj_pt_cut) {
	sjet3[nsub][1] = j->px();
	sjet3[nsub][2] = j->py();
	sjet3[nsub][3] = j->pz();
	sjet3[nsub][0] = j->e();
	vector<double> ys = jets.getYSubJet(*j);
	if (ys.size()>3) {
	  ysubsj[nsub][0] = ys.at(0);
	  ysubsj[nsub][1] = ys.at(1);
	  ysubsj[nsub][2] = ys.at(2);
	  ysubsj[nsub][3] = ys.at(3);
	} else {
	  ysubsj[nsub][0] = 0;
	  ysubsj[nsub][1] = 0;
	  ysubsj[nsub][2] = 0;
	  ysubsj[nsub][3] = 0;
	}
	nsub++;
      }
      njet++;
    }
  }

  // Loop over leptons
  nlep=0;
  for (ParticleVector::const_iterator p = cl.chargedLeptons().begin(); p != cl.chargedLeptons().end(); ++p) {
    if (p->getMomentum().perp()>_lepton_pt_cut) {
      vlep[nlep][1] = p->getMomentum().px();
      vlep[nlep][2] = p->getMomentum().py();
      vlep[nlep][3] = p->getMomentum().pz();
      vlep[nlep][0] = p->getMomentum().e();
      nlep++;
    }
  }

  // Total/missing energy.  
  esumr[1] = tvm.getMomentum().px();
  esumr[2] = tvm.getMomentum().py();
  esumr[3] = tvm.getMomentum().pz();
  esumr[0] = tvm.getMomentum().e();

  // Finished...
  log << Log::DEBUG << "Finished analyzing" << endl;

  rivetTree->Fill();
#endif
}


// Finalize
void ExampleTree::finalize() { 
#ifdef HAVE_ROOT
  // Write the tree to file.
  rivetTree->Write();
#endif
}

