// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Beams.hh"
using namespace std;

namespace Rivet {


  bool validBeams(const ParticlePair& beams) {
    return beams.first.pid() != 0 && beams.second.pid() != 0;
  }


  ParticlePair beams(const Event& e) {

    // Return a null pair if there aren't enough particles for there to be beams
    if (HepMCUtils::particles_size(e.genEvent()) >= 2) {

      // First try the official way: ask the GenEvent for the beam pointers
      std::pair<ConstGenParticlePtr,ConstGenParticlePtr> thebeams = HepMCUtils::beams(e.genEvent());
      if ( thebeams.first && thebeams.second ) {
        return ParticlePair{thebeams.first, thebeams.second};
      }
      // Ok, that failed: let's find the status = 4 particles by hand
      const vector<Particle> pstat4s = e.allParticles([](const Particle& p){ return p.genParticle()->status() == 4; });
      if (pstat4s.size() >= 2) {
        return ParticlePair{pstat4s[0], pstat4s[1]};
      }

      // Hmm, this sucks. Last guess is that barcodes 1 and 2 are the beams
      /// @todo There are no barcodes in HepMC3. Fall back to guessing from particle.id() rather than barcode?
      #ifndef RIVET_ENABLE_HEPMC_3
      if (e.genEvent()->barcode_to_particle(1) && e.genEvent()->barcode_to_particle(2)) {
        return ParticlePair{e.genEvent()->barcode_to_particle(1), e.genEvent()->barcode_to_particle(2)};
      }
      #endif
    }
    // else cerr << "EVENT TOO SMALL TO HAVE BEAMS!!!" << endl;

    // Give up: return null beams
    return ParticlePair(); //{Particle(), Particle()}; //< beam PIDs = ANY,ANY
  }


  double sqrtS(const FourMomentum& pa, const FourMomentum& pb) {
    const double mom1 = pa.pz();
    const double e1 = pa.E();
    const double mom2 = pb.pz();
    const double e2 = pb.E();
    const double sqrts = sqrt( sqr(e1+e2) - sqr(mom1+mom2) );
    return sqrts;
  }

  double sqrtS(double ea, double eb) {
    const double s = sqr(ea+eb) - sqr(ea-eb);
    const double sqrts = sqrt(s);
    return sqrts;
  }

  double asqrtS(const FourMomentum& pa, const FourMomentum& pb) {
    const static double MNUCLEON = 939*MeV; //< nominal nucleon mass
    return sqrtS(pa/(pa.mass()/MNUCLEON), pb/(pb.mass()/MNUCLEON));
  }

  double asqrtS(const ParticlePair& beams) {
    return sqrtS(beams.first.mom()/nuclA(beams.first), beams.second.mom()/nuclA(beams.second));
  }


  FourMomentum acmsBoostVec(const FourMomentum& pa, const FourMomentum& pb) {
    const static double MNUCLEON = 939*MeV; //< nominal nucleon mass
    const double na = pa.mass()/MNUCLEON, nb = pb.mass()/MNUCLEON;
    return cmsBoostVec(pa/na, pb/nb);
  }

  FourMomentum acmsBoostVec(const ParticlePair& beams) {
    return cmsBoostVec(beams.first.mom()/nuclA(beams.first), beams.second.mom()/nuclA(beams.second));
  }


  Vector3 cmsBetaVec(const FourMomentum& pa, const FourMomentum& pb) {
    // const Vector3 rtn = (pa.p3() + pb.p3()) / (pa.E() + pb.E());
    const Vector3 rtn = (pa + pb).betaVec();
    return rtn;
  }

  Vector3 acmsBetaVec(const FourMomentum& pa, const FourMomentum& pb) {
    const static double MNUCLEON = 939*MeV; //< nominal nucleon mass
    const Vector3 rtn = cmsBetaVec(pa/(pa.mass()/MNUCLEON), pb/(pb.mass()/MNUCLEON));
    return rtn;
  }

  Vector3 acmsBetaVec(const ParticlePair& beams) {
    const Vector3 rtn = cmsBetaVec(beams.first.mom()/nuclA(beams.first), beams.second.mom()/nuclA(beams.second));
    return rtn;
  }


  Vector3 cmsGammaVec(const FourMomentum& pa, const FourMomentum& pb) {
    // const Vector3 rtn = (pa + pb).gammaVec();
    const double gamma = (pa.E() + pb.E()) / sqrt( sqr(pa.mass()) + sqr(pb.mass()) + 2*(pa.E()*pb.E() - dot(pa.p3(), pb.p3())) );
    const Vector3 rtn = gamma * (pa.p3() + pb.p3()).unit();
    return rtn;
  }

  Vector3 acmsGammaVec(const FourMomentum& pa, const FourMomentum& pb) {
    const static double MNUCLEON = 939*MeV; //< nominal nucleon mass
    Vector3 rtn = cmsGammaVec(pa/(pa.mass()/MNUCLEON), pb/(pb.mass()/MNUCLEON));
    return rtn;
  }

  Vector3 acmsGammaVec(const ParticlePair& beams) {
    Vector3 rtn = cmsGammaVec(beams.first.mom()/nuclA(beams.first), beams.second.mom()/nuclA(beams.second));
    return rtn;
  }


  LorentzTransform cmsTransform(const FourMomentum& pa, const FourMomentum& pb) {
    /// @todo Automatically choose to construct from beta or gamma according to which is more precise?
    return LorentzTransform::mkFrameTransformFromGamma(cmsGammaVec(pa, pb));
  }

  LorentzTransform acmsTransform(const FourMomentum& pa, const FourMomentum& pb) {
    /// @todo Automatically choose to construct from beta or gamma according to which is more precise?
    return LorentzTransform::mkFrameTransformFromGamma(acmsGammaVec(pa, pb));
  }

  LorentzTransform acmsTransform(const ParticlePair& beams) {
    return LorentzTransform::mkFrameTransformFromGamma(acmsGammaVec(beams));
  }


}
