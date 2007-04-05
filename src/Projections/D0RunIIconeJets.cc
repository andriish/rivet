// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/D0RunIIconeJets.hh"
#include "Rivet/Projections/Cmp.hh"
#include "Rivet/RivetCLHEP.hh"

#include "Rivet/Tools/D0RunIIcone/inline_maths.h"
#include "Rivet/Tools/D0RunIIcone/HepEntity.h"

using namespace Rivet;

/* The compare method is incomplete? */

int D0RunIIconeJets::compare(const Projection & p) const {
  const D0RunIIconeJets & other = dynamic_cast<const D0RunIIconeJets &>(p);
  return pcmp(*fsproj, *other.fsproj) || cmp(cone_radius, other.cone_radius) ||
    cmp(min_jet_Et, other.min_jet_Et) || cmp(split_ratio, other.split_ratio);
}

void D0RunIIconeJets::project(const Event & e) {
  //vector<KtJet::KtLorentzVector> vecs;
  particlelist = new list<const HepEntity*>; //privately declared

  // Project into final state
  const FinalState& fs = e.applyProjection(*fsproj);
  
  //double SET=0.;
  // Store 4 vector data about each particle into list
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    LorentzVector fv = p->getMomentum();
    // store the FourVector in the KtLorentzVector form
    //KtJet::KtLorentzVector ktlv(fv.px(), fv.py(), fv.pz(), fv.e());
    const HepEntity* listelement = new HepEntity(fv.e(), fv.px(), fv.py(), fv.pz());
    particlelist->push_back(listelement);
    //SET += sqrt(sqr(fv.px())+sqr(fv.py()));
}

  //cout << "Final State Particles (" << fs.particles().size() << ") SET = " << SET << endl;

  //pktev_ = new KtJet::KtEvent(vecs, type_, angle_, recom_, rparameter_);
  float Item_ET_Threshold = 0.;
  float* Item_ET_Threshold_ptr = &Item_ET_Threshold;

  //jets = pointer to list of type HepEntity 
  if (jets->size()>0) jets->clear(); //be sure to have no jets from old event(s)
  algo->makeClusters(*jets, *particlelist, Item_ET_Threshold); //turn the crank!!!
 
  /*  
  list<HepEntity>::iterator it;
  cout << "Number of jets = " << jets.size() << endl;
  for (it=jets->begin(); it!=jets->end(); ++it) {
    cout << "jet: E=" << (*it).E << " pT=" << (*it).pT() << " y=" << (*it).y() << " phi=" << (*it).phi() << endl; 
  }
  */

}

RivetInfo D0RunIIconeJets::getInfo() const {
  return Projection::getInfo() + fsproj->getInfo();
}


