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


  double asqrtS(const FourMomentum& pa, const FourMomentum& pb) {
    const static double MNUCLEON = 939*MeV; //< nominal nucleon mass
    return sqrtS(pa/(pa.mass()/MNUCLEON), pb/(pb.mass()/MNUCLEON));
  }

  double asqrtS(const ParticlePair& beams) {
    return sqrtS(beams.first.mom()/nuclA(beams.first), beams.second.mom()/nuclA(beams.second));
  }


  Vector3 cmsBoostBeta(const FourMomentum& pa, const FourMomentum& pb) {
    Vector3 rtn = (pa.p3() + pb.p3()) / (pa.E() + pb.E()); ///< @todo Is the vector part of this correct in general?
    return rtn;
  }

  Vector3 acmsBoostBeta(const FourMomentum& pa, const FourMomentum& pb) {
    const static double MNUCLEON = 939*MeV; //< nominal nucleon mass
    Vector3 rtn = cmsBoostBeta(pa/(pa.mass()/MNUCLEON), pb/(pb.mass()/MNUCLEON));
    return rtn;
  }

  Vector3 acmsBoostBeta(const ParticlePair& beams) {
    Vector3 rtn = cmsBoostBeta(beams.first.mom()/nuclA(beams.first), beams.second.mom()/nuclA(beams.second));
    return rtn;
  }


  Vector3 cmsBoostGamma(const FourMomentum& pa, const FourMomentum& pb) {
    const double gamma = ( sqr(pa.mass()) + sqr(pb.mass()) + 2*(pa.E()*pb.E() - dot(pa.p3(), pb.p3())) ) / sqr(pa.E() + pb.E());
    Vector3 rtn = gamma * cmsBoostBeta(pa, pb).unit(); ///< @todo Duh. There must be a better way to get the unit vector of boost direction...
    return rtn;
  }

  Vector3 acmsBoostGamma(const FourMomentum& pa, const FourMomentum& pb) {
    const static double MNUCLEON = 939*MeV; //< nominal nucleon mass
    Vector3 rtn = cmsBoostGamma(pa/(pa.mass()/MNUCLEON), pb/(pb.mass()/MNUCLEON));
    return rtn;
  }

  Vector3 acmsBoostGamma(const ParticlePair& beams) {
    Vector3 rtn = cmsBoostGamma(beams.first.mom()/nuclA(beams.first), beams.second.mom()/nuclA(beams.second));
    return rtn;
  }


  LorentzTransform cmsTransform(const FourMomentum& pa, const FourMomentum& pb) {
    /// @todo Automatically choose to construct from beta or gamma according to which is more precise?
    return LorentzTransform::mkFrameTransformFromGamma(cmsBoostGamma(pa, pb));
  }

  LorentzTransform acmsTransform(const FourMomentum& pa, const FourMomentum& pb) {
    /// @todo Automatically choose to construct from beta or gamma according to which is more precise?
    return LorentzTransform::mkFrameTransformFromGamma(acmsBoostGamma(pa, pb));
  }

  LorentzTransform acmsTransform(const ParticlePair& beams) {
    return LorentzTransform::mkFrameTransformFromGamma(acmsBoostGamma(beams));
  }



  /////////////////////////////////////////////



  void Beam::project(const Event& e) {
    _theBeams = Rivet::beams(e);
    MSG_DEBUG("Beam particles = " << _theBeams << " => sqrt(s) = " << sqrtS()/GeV << " GeV");
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
