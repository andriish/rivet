// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/HepEx0112029.hh"
using namespace Rivet;

#include "Rivet/RivetAIDA.hh"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////


// Book histograms
void HepEx0112029::init() {
  /// @todo This doesn't seem to correspond to the plots in the paper (SPIRES 4730372)
  _histJetEt1 = bookHistogram1D("JetET1", "Jet transverse energy", 11, 14.0, 75.0);
}


// Do the analysis
void HepEx0112029::analyze(const Event& event) {
  Log& log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;

  // Analyse and print some info
  const KtJets& jets = event.applyProjection(_ktjetsproj);
  const int nj = jets.getNJets();
  log << Log::INFO << "Jet multiplicity = " << nj << endl;

  // Fill histograms
  vector<KtJet::KtLorentzVector> jetList = jets.getJetsEt();
  for (vector<KtJet::KtLorentzVector>::iterator j = jetList.begin(); j != jetList.end(); ++j) {
    _histJetEt1->fill(j->perp(), event.weight() );
  }
  
  // Finished
  log << Log::DEBUG << "Finished analyzing" << endl;
}


// Finalize
void HepEx0112029::finalize() 
{ }
