// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FinalStateHCM class.
//

#include "Rivet/Projections/FinalStateHCM.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;

FinalStateHCM::~FinalStateHCM() {}

int FinalStateHCM::compare(const Projection & p) const {
  const FinalStateHCM & other =
    dynamic_cast<const FinalStateHCM &>(p);
  return pcmp(lepton, other.lepton) ||
    pcmp(kinematics, other.kinematics) || pcmp(fsproj, other.fsproj);
}

void FinalStateHCM::project(const Event & e) {
  const DISLepton & dislep = e(lepton);
  const DISKinematics & diskin = e(kinematics);
  const FinalState & fs = e(fsproj);
  theParticles.clear();
  theParticles.reserve(fs.particles().size());
  for ( int i = 0, N = fs.particles().size(); i < N; ++i )
    if ( fs.particles()[i].original != dislep.out().original ) {
      theParticles.push_back(fs.particles()[i]);
      theParticles[i].momentum *= diskin.boostHCM();
    }
}

RivetInfo FinalStateHCM::getInfo() const {
  return Projection::getInfo() + lepton.getInfo() +
    kinematics.getInfo() + fsproj.getInfo();
}

