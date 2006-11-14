// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CentralEtHCM class.
//

#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;
using namespace std;

CentralEtHCM::~CentralEtHCM() {}

int CentralEtHCM::compare(const Projection & p) const {
  const CentralEtHCM & other =
    dynamic_cast<const CentralEtHCM &>(p);
  return pcmp(fshcm, other.fshcm);
}

void CentralEtHCM::project(const Event & e) {
  const FinalStateHCM & fs = e(fshcm);
  sumet = 0.0;
  for ( int i = 0, N = fs.particles().size(); i < N; ++i ) {
    if ( abs(fs.particles()[i].momentum.rapidity()) < 0.5 )
      sumet += fs.particles()[i].momentum.et();
  }
}

RivetInfo CentralEtHCM::getInfo() const {
  return Projection::getInfo() + fshcm.getInfo();
}

