// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {

  int DISLepton::compare(const Projection& p) const {
    const DISLepton& other = pcast<DISLepton>(p);
    return
      mkNamedPCmp(other, "Beam") || 
      mkNamedPCmp(other, "FS") || 
      cmp(_idin, other._idin) || 
      cmp(_idout, other._idout);
  }


  void DISLepton::project(const Event& e) {
    const ParticlePair& inc = applyProjection<Beam>(e, "Beam").getBeams();
    const FinalState& fs = applyProjection<FinalState>(e, "FS");

    bool allowAnti = (_idin * _idout < 0);
    if ( _idin == inc.first.getPdgId() || (allowAnti && _idin == -inc.first.getPdgId()) ) {
      _incoming = inc.first;
    } else if ( _idin == inc.second.getPdgId() || (allowAnti && _idin == -inc.second.getPdgId()) ) {
      _incoming = inc.second;
    } else {
      throw	Error("DISLepton projector could not find the correct beam. ");
    }

    double emax = 0.0;
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      if ( ( _idout == p->getPdgId() || (allowAnti && _idout == -p->getPdgId()) ) && 
           p->momentum().E() > emax ) {
        /// @todo change this to a correct way of finding the scattered lepton.
        emax = p->momentum().E();
        _outgoing = *p;
      }
    }
    
    if (!_outgoing.hasHepMCParticle()) {
      throw Error("DISLepton projector could not find the scattered lepton.");
    }
  }
  
}

