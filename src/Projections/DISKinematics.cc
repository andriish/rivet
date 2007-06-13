// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/Cmp.hh"


namespace Rivet {

  void DISKinematics::project(const Event& e) {
    const DISLepton& dislep = e.applyProjection(*_lepton);
    const ParticlePair& inc = e.applyProjection(*_beams).getBeams();
    Particle hadron;

    if ( inc.second.getPdgId() == _idhad ) hadron = inc.second;
    else if ( inc.first.getPdgId() == _idhad ) hadron = inc.first;
    else throw runtime_error("DISKinematics projector could not find "
			     "the correct beam.");

    if ( &(dislep.in().getHepMCParticle()) ==
	 &(hadron.getHepMCParticle()) )
      throw runtime_error("DISKinematics projector could not find the "
			  "correct beam.");

    LorentzVector pgam(dislep.in().getMomentum() - dislep.out().getMomentum());
    _theQ2 = -pgam.m2();
    LorentzVector tothad(pgam + hadron.getMomentum());
    _theW2 = tothad.m2();
    _theX = Q2()/(2.0*pgam*hadron.getMomentum());
    _hcm.set(-tothad.boostVector());
    pgam *= _hcm;
    _hcm.rotateX(-pgam.theta());
    LorentzVector plep = dislep.out().getMomentum();
    plep *= _hcm;
    _hcm.rotateZ(-plep.phi());
    _breit = _hcm;
    double bz = (pgam.rho() - 0.5*x()*sqrt(W2()))/(pgam.e() + 0.5*x()*sqrt(W2()));
    _breit.boostZ(-bz);
    // @todo *** ATTENTION *** These rotations should be checked.
  }


  int DISKinematics::compare(const Projection & p) const {
    const DISKinematics& other = dynamic_cast<const DISKinematics&>(p);
    return \
      pcmp(*_lepton, *other._lepton) || 
      pcmp(*_beams, *other._beams) || 
      cmp(_idhad, other._idhad);
  }


}
