// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISKinematics class.
//

#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;

DISKinematics::~DISKinematics() {}

void DISKinematics::project(const Event& e) {
  const DISLepton & dislep = e.applyProjection(lepton);
  const ParticlePair & inc = e.applyProjection(beams)();
  Particle hadron;
  if ( inc.second.id == idhad ) hadron = inc.second;
  else
    throw
      runtime_error("DISKinematics projector could not find the correct beam.");
  // *** ATTENTION *** Should we have our own exception classes?
    
  if ( dislep.in().original == hadron.original )
    throw
      runtime_error("DISKinematics projector could not find the correct beam.");
  // *** ATTENTION *** Should we have our own exception classes?

  LorentzVector pgam(dislep.in().momentum - dislep.out().momentum);
  theQ2 = -pgam.m2();
  LorentzVector tothad(pgam + hadron.momentum);
  theW2 = tothad.m2();
  theX = Q2()/(2.0*pgam*hadron.momentum);
  hcm.set(-tothad.boostVector());
  pgam *= hcm;
  hcm.rotateX(-pgam.theta());
  LorentzVector plep = dislep.out().momentum;
  plep *= hcm;
  hcm.rotateZ(-plep.phi());
  breit = hcm;
  double bz = (pgam.rho() - 0.5*x()*sqrt(W2()))/(pgam.e() + 0.5*x()*sqrt(W2()));
  breit.boostZ(-bz);
  // *** ATTENTION *** These rotations should be checked.
}

int DISKinematics::compare(const Projection & p) const {
  const DISKinematics & other = dynamic_cast<const DISKinematics &>(p);
  return pcmp(lepton, other.lepton) || pcmp(beams, other.beams) ||
    cmp(idhad, other.idhad);
}

RivetInfo DISKinematics::getInfo() const {
  return Projection::getInfo() + beams.getInfo() + lepton.getInfo();
}

