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
  _lepton_pt_cut = 20;

  treeFileName = "rivetTree.root";

  // Create a file for the Tree
  treeFile = new TFile(treeFileName,"recreate");

  // book the ntuple.
  rivetTree = new TTree("Rivet Tree","Rivet Example Tree");

  // event number 
  rivetTree->Branch("nevt",&nevt,"nevt/I");

  // 
  rivetTree->Branch("nw",&nw,"nw/I");
  rivetTree->Branch("wtype",&wtype,"wtype[nw][2]/I");
  rivetTree->Branch("wvec",&wvec,"wtype[nw][4]/F");
  rivetTree->Branch("hvec",&hvec,"hvec[4]/F");

  rivetTree->Branch("njet",&njet,"njet/I");
  rivetTree->Branch("ptjet",&ptjet,"ptjet[njet]/F");
  rivetTree->Branch("etajet",&etajet,"etajet[njet]/F");
  rivetTree->Branch("phijet",&phijet,"phijet[njet]/F");
  rivetTree->Branch("vjet",&vjet,"vjet[njet][4]/F");

  rivetTree->Branch("nsub",&nsub,"nsub/I");
  rivetTree->Branch("sjet3",&sjet3,"sjet3[nsub][4]/F");
  rivetTree->Branch("ysubsj",&ysubsj,"ysubsj[nsub]/F");

  rivetTree->Branch("tjet",&tjet,"tjet[2][4]/F");

  rivetTree->Branch("nlep",&nlep,"nlep/I");
  rivetTree->Branch("ptlep",&ptlep,"ptlep[nlep]/F");
  rivetTree->Branch("etalep",&etalep,"etalep[nlep]/F");
  rivetTree->Branch("philep",&philep,"philep[nlep]/F");
  rivetTree->Branch("vlep",&vlep,"vlep[nlep][4]/F");

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

  // @todo W bosons in the event
  nw = 0;

  // Jets.
  const KtJets& jets = event.applyProjection(p_ktjets);

  // Leptons
  const ChargedLeptons& cl = event.applyProjection(p_chargedleptons);

  // Missing Et/total energy
  const TotalVisibleMomentum& tvm = event.applyProjection(*p_totalvisiblemomentum);

  // Vector bosons.


  // Get the jets in decreasing ET order.
  vector<KtJet::KtLorentzVector> jetList = jets.getJetsEt();
  njet = 0;
  for (vector<KtJet::KtLorentzVector>::iterator j = jetList.begin(); j != jetList.end(); ++j) {
    if (j->perp()>_jet_pt_cut) {
      ptjet[njet] = j->perp();
      etajet[njet] = j->eta();
      phijet[njet] = j->phi();
      vjet[njet][0] = j->px();
      vjet[njet][1] = j->py();
      vjet[njet][2] = j->pz();
      vjet[njet][3] = j->e();
      njet++;
    }
  }

  // Loop over leptons
  nlep=0;
  for (ParticleVector::const_iterator p = cl.chargedLeptons().begin(); p != cl.chargedLeptons().end(); ++p) {
    if (p->getMomentum().perp()>_lepton_pt_cut) {
      ptlep[nlep] = p->getMomentum().perp();
      etalep[nlep] = p->getMomentum().eta();
      philep[nlep] = p->getMomentum().phi();
      vlep[nlep][0] = p->getMomentum().px();
      vlep[nlep][1] = p->getMomentum().py();
      vlep[nlep][2] = p->getMomentum().pz();
      vlep[nlep][3] = p->getMomentum().e();
      nlep++;
    }
  }

  // Total/missing energy.  
  esumr[0] = tvm.getMomentum().px();
  esumr[1] = tvm.getMomentum().py();
  esumr[2] = tvm.getMomentum().pz();
  esumr[3] = tvm.getMomentum().e();

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

