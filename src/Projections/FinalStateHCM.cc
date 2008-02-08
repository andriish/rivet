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
    const GenParticle* dislepGP = &( dislep.out().getHepMCParticle() );

    const DISKinematics& diskin = e.applyProjection(_kinematics);
    const LorentzTransform hcmboost = diskin.boostHCM();

    const FinalState& fs = e.applyProjection(_fsproj);

    // Fill the particle list with all particles _other_ than the DIS scattered
    // lepton, with momenta boosted into the HCM frame.
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());
    for (size_t i = 0, N = fs.particles().size(); i < N; ++i) {
      const GenParticle* loopGP = &( fs.particles()[i].getHepMCParticle() );
      /// @todo Maybe Particle should have equiv() and same() methods?
      if (loopGP != dislepGP) {
        Particle tmpP = fs.particles()[i];
        const FourMomentum hcmMom = hcmboost.transform(tmpP.getMomentum());
        tmpP.setMomentum(hcmMom);
        _theParticles.push_back(tmpP);
      }
    }
  }

}
