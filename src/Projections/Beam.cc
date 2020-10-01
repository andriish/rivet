// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  void Beam::project(const Event& e) {
    _theBeams = Rivet::beams(e);
    MSG_DEBUG("Beam particles = " << _theBeams << " => sqrt(s) = " << sqrtS()/GeV << " GeV");
  }


  FourVector Beam::pv() const {
    RivetHepMC::FourVector v1, v2;
    const ParticlePair bpair = beams();
    if (bpair.first.genParticle() && bpair.first.genParticle()->end_vertex())
      v1 = bpair.first.genParticle()->end_vertex()->position();
    if (bpair.second.genParticle() && bpair.second.genParticle()->end_vertex())
      v2 = bpair.second.genParticle()->end_vertex()->position();
    const FourVector rtn = (v1 == v2) ? FourVector(v1.t(), v1.x(), v1.y(), v1.z()) : FourVector();
    MSG_DEBUG("Beam PV 4-position = " << rtn);
    return rtn;
  }


}
