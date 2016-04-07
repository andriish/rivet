// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  ParticlePair beams(const Event& e) {
    assert(e.genEvent()->particles_size() >= 2);
    if (e.genEvent()->valid_beam_particles()) {
      pair<HepMC::GenParticle*, HepMC::GenParticle*> beams = e.genEvent()->beam_particles();
      assert(beams.first && beams.second);
      return ParticlePair{beams.first, beams.second};
    } else if (e.genEvent()->barcode_to_particle(1) && e.genEvent()->barcode_to_particle(2)) {
      return ParticlePair{e.genEvent()->barcode_to_particle(1), e.genEvent()->barcode_to_particle(2)};
    }
    return ParticlePair{Particle(PID::ANY, FourMomentum()), Particle(PID::ANY, FourMomentum())};
  }


  double sqrtS(const FourMomentum& pa, const FourMomentum& pb) {
    const double mom1 = pa.pz();
    const double e1 = pa.E();
    const double mom2 = pb.pz();
    const double e2 = pb.E();
    double sqrts = sqrt( sqr(e1+e2) - sqr(mom1+mom2) );
    return sqrts;
  }


  double asqrtS(const ParticlePair& beams) {
    return sqrtS(beams) / (nuclA(beams.first) + nuclA(beams.second));
  }



  /////////////////////////////////////////////



  void Beam::project(const Event& e) {
    _theBeams = Rivet::beams(e);
    MSG_DEBUG("Beam particles = " << _theBeams << " => sqrt(s) = " << sqrtS()/GeV << " GeV");
  }


  double Beam::sqrtS() const {
    return Rivet::sqrtS(beams());
  }


  FourVector Beam::pv() const {
    HepMC::FourVector v1, v2;
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
