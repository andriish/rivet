// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "KtJet/KtJet.h"
#include "HepMC/GenEvent.h"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/KtJets.hh"
#include <fstream>

using namespace Rivet;

KtJets::~KtJets() {}

/* The compare method is incomplete? */

int KtJets::compare(const Projection & p) const {
  //const KtJets & other = dynamic_cast<const KtJets &>(p);
  return 0;
}

void KtJets::getJetsE() {   
  for(std::vector<KtJet::KtLorentzVector>::iterator vecIt = jetVecs_.begin(); vecIt!= jetVecs_.begin()+4; ++vecIt){
    //ofstream book_file("BookInfo.cc");
    std::cout<<vecIt->px()<<", "<<vecIt->py()<<", "<<vecIt->pz()<<", "<<vecIt->e()<<endl; 
    //book_file<<vecIt->px()<<", "<<vecIt->py()<<", "<<vecIt->pz()<<", "<<vecIt->e()<<endl;
  }
}

void KtJets::project(const Event & e) {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);
  
  std::vector<KtJet::KtLorentzVector> vecs;
  
  // Clear counters
  NJets_ = 0;
    
  // Project into final state
  const FinalState& fs = e.addProjection(fsproj);
  
  // Store 4 vector data about each particle into vecs
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    HepMC::FourVector fv2 = p->momentum;
    // store the FourVector in the KtLorentzVector form
    KtJet::KtLorentzVector ktlv2(fv2.px(), fv2.py(), fv2.pz(), fv2.e());
    vecs.push_back(ktlv2);
  }
  
  KtJet::KtEvent ktev(vecs, 4, 3, 1, 1.0);
  NJets_ = ktev.getNJets();
  jetVecs_ = ktev.getJetsE();
}
