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
    const ParticlePair& inc = applyProjection<Beam>(e, "Beam").beams();
    const bool allowAnti = (_idin * _idout < 0);
    if ( _idin == inc.first.pdgId() || (allowAnti && _idin == -inc.first.pdgId()) ) {
      _incoming = inc.first;
    } else if ( _idin == inc.second.pdgId() || (allowAnti && _idin == -inc.second.pdgId()) ) {
      _incoming = inc.second;
    } else {
      throw	Error("DISLepton projector could not find the correct beam. ");
    }

    double emax = 0.0;
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    foreach (const Particle& p, fs.particles()) {
      if ( ( _idout == p.pdgId() || (allowAnti && _idout == -p.pdgId()) ) && 
           p.momentum().E() > emax ) {
        /// @todo change this to a correct way of finding the scattered lepton.
        emax = p.momentum().E();
        _outgoing = p;
      }
    }
    
    if (!_outgoing.hasGenParticle()) {
      throw Error("DISLepton projector could not find the scattered lepton.");
    }
  }

  
}
