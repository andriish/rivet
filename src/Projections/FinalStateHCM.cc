// -*- C++ -*-

#include "Rivet/Projections/FinalStateHCM.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;


int FinalStateHCM::compare(const Projection & p) const {
  const FinalStateHCM & other =
    dynamic_cast<const FinalStateHCM &>(p);
  return pcmp(*lepton, *other.lepton) ||
    pcmp(*kinematics, *other.kinematics) || pcmp(*fsproj, *other.fsproj);
}

void FinalStateHCM::project(const Event& e) {
  const DISLepton & dislep = e.applyProjection(*lepton);
  const DISKinematics & diskin = e.applyProjection(*kinematics);
  const FinalState & fs = e.applyProjection(*fsproj);
  theParticles.clear();
  theParticles.reserve(fs.particles().size());
  for ( int i = 0, N = fs.particles().size(); i < N; ++i ) {
    if ( fs.particles()[i].getHepMCParticle() != dislep.out().getHepMCParticle() ) {
      theParticles.push_back(fs.particles()[i]);
      theParticles[i].getMomentum() *= diskin.boostHCM();
    }
  }
}
