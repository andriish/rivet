// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/Cmp.hh"


namespace Rivet {

  void DISLepton::project(const Event& e) {
    const ParticlePair& inc = e.applyProjection(*_beams)();
      // @todo Should this be _idin == -_idout?
    if ( (_idin == -_idin && abs(inc.first.getPdgId()) == abs(_idin)) || inc.first.getPdgId() ) {
      _incoming = inc.first;
    } else {
      throw runtime_error("DISLepton projector could not find the correct beam. ");
    }

    double emax = 0.0;
    for (GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
          pi != e.genEvent().particles_end(); ++pi) {
      // @todo Should this be _idin == -_idout?
      if ((_idin == -_idin && abs((*pi)->pdg_id()) == abs(_idout)) ||
          (*pi)->pdg_id() == _idout && (*pi)->momentum().e() > emax) {
        /// @todo This is probably not the correct way to select the scattered lepton
        emax = (*pi)->momentum().e();
        _outgoing = Particle(**pi);
      }
    }
  }


  int DISLepton::compare(const Projection& p) const {
    const DISLepton& other = dynamic_cast<const DISLepton&>(p);
    return \
      pcmp(*_beams, *other._beams) || 
      cmp(_idin, other._idin) || 
      cmp(_idout, other._idout);
  }

}

