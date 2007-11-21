// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {

  void DISKinematics::project(const Event& e) {
    const DISLepton& dislep = e.applyProjection(_lepton);
    const ParticlePair& inc = e.applyProjection(_beams).getBeams();
    Particle hadron;

    if ( inc.second.getPdgId() == _idhad ) hadron = inc.second;
    else if ( inc.first.getPdgId() == _idhad ) hadron = inc.first;
    else throw runtime_error("DISKinematics projector could not find the correct beam.");

    if ( &(dislep.in().getHepMCParticle()) == &(hadron.getHepMCParticle()) ) {
      throw runtime_error("DISKinematics projector could not find the correct beam.");
    }

    const FourMomentum phad = hadron.getMomentum();
    FourMomentum pgam(dislep.in().getMomentum() - dislep.out().getMomentum());
    _theQ2 = -pgam.mass2();
    FourMomentum tothad(pgam + phad);
    _theW2 = tothad.mass2();
    _theX = Q2()/(2 * pgam * phad);
    _hcm.setBoost(-boostVector(tothad));
    pgam = _hcm.transform(pgam);
    _hcm.rotate(Vector3::X(), -pgam.polarAngle());
    FourMomentum plep = dislep.out().getMomentum();
    plep = _hcm.transform(plep);
    _hcm.rotate(Vector3::Z(), -plep.azimuthalAngle());
    _breit = _hcm;
    const double bz = ( pgam.pT() - 0.5*x()*sqrt(W2()) ) / ( pgam.E() + 0.5*x()*sqrt(W2()) );
    _breit.setBoost(Vector3::Z() * -bz);
    /// @todo *** ATTENTION *** These rotations should be checked.
  }


  int DISKinematics::compare(const Projection & p) const {
    const DISKinematics& other = dynamic_cast<const DISKinematics&>(p);
    return \
      pcmp(_lepton, other._lepton) || 
      pcmp(_beams, other._beams) || 
      cmp(_idhad, other._idhad);
  }


}
