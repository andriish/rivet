// -*- C++ -*-

#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/Projections/Cmp.hh"
#include "Rivet/RivetCLHEP.hh"

using namespace Rivet;
using namespace std;


int CentralEtHCM::compare(const Projection & p) const {
  const CentralEtHCM & other =
    dynamic_cast<const CentralEtHCM &>(p);
  return pcmp(*fshcm, *other.fshcm);
}

void CentralEtHCM::project(const Event & e) {
  const FinalStateHCM & fs = e.applyProjection(*fshcm);
  sumet = 0.0;
  for ( int i = 0, N = fs.particles().size(); i < N; ++i ) {
    if ( abs(fs.particles()[i].getMomentum().rapidity()) < 0.5 )
      sumet += fs.particles()[i].getMomentum().et();
  }
}

// RivetInfo CentralEtHCM::getInfo() const {
//   return Projection::getInfo() + fshcm->getInfo();
// }

