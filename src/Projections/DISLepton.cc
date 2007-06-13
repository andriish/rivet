// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/RivetCLHEP.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/Cmp.hh"


namespace Rivet {

  void DISLepton::project(const Event& e) {
    const ParticlePair& inc = e.applyProjection(*_beams)();
    const FinalState& fs = e.applyProjection(*_fsproj);

    bool allowanti = ( _idin*_idout < 0 );

    if ( _idin == inc.first.getPdgId() ||
	 ( allowanti && _idin == -inc.first.getPdgId() ) )
      _incoming = inc.first;
    else if ( _idin == inc.second.getPdgId() ||
	 ( allowanti && _idin == -inc.second.getPdgId() ) )
      _incoming = inc.second;
    else
      throw
	runtime_error("DISLepton projector could not find the correct beam. ");

    double emax = 0.0;
    for (ParticleVector::const_iterator p = fs.particles().begin();
	 p != fs.particles().end(); ++p) {
      if ( ( _idout == p->getPdgId() ||
	     allowanti && _idout == -p->getPdgId() ) &&
	   p->getMomentum().e() > emax ) {
	// @todo change this to a correct way of finding the scattered lepton.
	emax = p->getMomentum().e();
	_outgoing = *p;
      }
    }

    if ( !_outgoing.hasHepMCParticle() )
      throw runtime_error("DISLepton projector could not find "
			  "the scattered lepton. ");


  }


  int DISLepton::compare(const Projection& p) const {
    const DISLepton& other = dynamic_cast<const DISLepton&>(p);
    return
      pcmp(*_beams, *other._beams) || 
      pcmp(*_fsproj, *other._fsproj) || 
      cmp(_idin, other._idin) || 
      cmp(_idout, other._idout);
  }

}

