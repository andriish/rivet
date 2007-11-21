// -*- C++ -*-

#include "Rivet/Projections/FinalStateHCM.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {

  int FinalStateHCM::compare(const Projection& p) const {
    const FinalStateHCM& other = dynamic_cast<const FinalStateHCM&>(p);
    return \
      pcmp(_lepton, other._lepton) ||
      pcmp(_kinematics, other._kinematics) || 
      pcmp(_fsproj, other._fsproj);
  }


  void FinalStateHCM::project(const Event& e) {
    const DISLepton& dislep = e.applyProjection(_lepton);
    const DISKinematics& diskin = e.applyProjection(_kinematics);
    const FinalState& fs = e.applyProjection(_fsproj);
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());
    for (size_t i=0, N=fs.particles().size(); i < N; ++i) {
      if ( &(fs.particles()[i].getHepMCParticle()) !=
	   &(dislep.out().getHepMCParticle()) ) {
        _theParticles.push_back(fs.particles()[i]);
        FourMomentum p_i = _theParticles[i].getMomentum();
        p_i = diskin.boostHCM().transform(p_i);
        _theParticles[i].setMomentum(p_i);
      }
    }
  }

}
