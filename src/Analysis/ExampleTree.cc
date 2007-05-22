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

  treeFileName = "rivetTree.root";

  // Create a file for the Tree
  treeFile = new TFile(treeFileName,"recreate");

  // book the ntuple.
  rivetTree = new TTree("Rivet Tree","Rivet Example Tree");
  rivetTree->Branch("nevt",&nevt,"nevt/I");

  rivetTree->Branch("nw",&nw,"nw/I");
  rivetTree->Branch("wtype",&wtype,"wtype[nw][2]/I");
  rivetTree->Branch("wvec",&wvec,"wtype[nw][4]/F");
  rivetTree->Branch("hvec",&hvec,"hvec[4]/F");

  rivetTree->Branch("njet",&njet,"njet/I");

  rivetTree->Branch("nsub",&nsub,"nsub/I");

#endif
}


// Do the analysis
void ExampleTree::analyze(const Event & event) {
#ifdef HAVE_ROOT
  Log log = getLog();
  log << Log::DEBUG << "Filling the ntuple" << endl;

  // Analyse and print some info
  const Multiplicity& m = event.applyProjection(p_mult);

  
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

