// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISKinematics class.
//

#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;


void DISKinematics::project(const Event& e) {
  const DISLepton & dislep = e.applyProjection(*lepton);
  const ParticlePair & inc = e.applyProjection(*beams)();
  Particle hadron;

  if (inc.second.getPdgId() == idhad) hadron = inc.second;
  else throw runtime_error("DISKinematics projector could not find the correct beam.");
  // @todo *** ATTENTION *** Should we have our own exception classes?
    
  if (dislep.in().getHepMCParticle() == hadron.getHepMCParticle())
    throw runtime_error("DISKinematics projector could not find the correct beam.");
  // @todo *** ATTENTION *** Should we have our own exception classes?

  LorentzVector pgam(dislep.in().getMomentum() - dislep.out().getMomentum());
  theQ2 = -pgam.m2();
  LorentzVector tothad(pgam + hadron.getMomentum());
  theW2 = tothad.m2();
  theX = Q2()/(2.0*pgam*hadron.getMomentum());
  hcm.set(-tothad.boostVector());
  pgam *= hcm;
  hcm.rotateX(-pgam.theta());
  LorentzVector plep = dislep.out().getMomentum();
  plep *= hcm;
  hcm.rotateZ(-plep.phi());
  breit = hcm;
  double bz = (pgam.rho() - 0.5*x()*sqrt(W2()))/(pgam.e() + 0.5*x()*sqrt(W2()));
  breit.boostZ(-bz);
  // @todo *** ATTENTION *** These rotations should be checked.
}

int DISKinematics::compare(const Projection & p) const {
  const DISKinematics & other = dynamic_cast<const DISKinematics &>(p);
  return pcmp(*lepton, *other.lepton) || pcmp(*beams, *other.beams) ||
    cmp(idhad, other.idhad);
}

// RivetInfo DISKinematics::getInfo() const {
//   return Projection::getInfo() + beams->getInfo() + lepton->getInfo();
// }
