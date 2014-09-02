// -*- C++ -*-
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @todo Not necessary in C++11, where constructors can be chained
  void TauFinder::_init(const FinalState& inputfs, DecayType decaytype) {
    setName("TauFinder");

    _taus.clear();

    addProjection(UnstableFinalState(), "UFS");

    //VetoedFinalState remainingFS;
    //remainingFS.addVetoOnThisFinalState(*this);
    //addProjection(remainingFS, "RFS");
  }


  void TauFinder::project(const Event& e) {
    clear();

    const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");
    foreach (const Particle& p, ufs.particles()) {
      if (p.abspid() != PID::TAU) continue;
      _taus.push_back(p);

    }
  }



}

