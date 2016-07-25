#ifndef RIVET_PARTICLEBASEUTILS_HH
#define RIVET_PARTICLEBASEUTILS_HH

#include "Rivet/ParticleBase.hh"

namespace Rivet {


  /// @name Particle classifying functors
  ///
  /// To be passed to any() or all() e.g. any(p.children(), HasPID(PID::MUON))
  //@{

  /// Base type for Particle -> bool functors
  struct BoolParticleBaseFunctor {
    virtual bool operator()(const ParticleBase& p) const = 0;
  };

  /// Transverse momentum greater-than functor
  struct PtGtr : public BoolParticleBaseFunctor {
    PtGtr(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() > ptcut; }
    double ptcut;
  };

  /// Transverse momentum less-than functor
  struct PtLess : public BoolParticleBaseFunctor {
    PtLess(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() < ptcut; }
    double ptcut;
  };


  /// Pseudorapidity greater-than functor
  struct EtaGtr : public BoolParticleBaseFunctor {
    EtaGtr(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() > etacut; }
    double etacut;
  };

  /// Pseudorapidity momentum less-than functor
  struct EtaLess : public BoolParticleBaseFunctor {
    EtaLess(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() < etacut; }
    double etacut;
  };

  /// Abs pseudorapidity greater-than functor
  struct AbsEtaGtr : public BoolParticleBaseFunctor {
    AbsEtaGtr(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() > absetacut; }
    double absetacut;
  };

  /// Abs pseudorapidity momentum less-than functor
  struct AbsEtaLess : public BoolParticleBaseFunctor {
    AbsEtaLess(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() < absetacut; }
    double absetacut;
  };


  /// Rapidity greater-than functor
  struct RapGtr : public BoolParticleBaseFunctor {
    RapGtr(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() > rapcut; }
    double rapcut;
  };

  /// Rapidity momentum less-than functor
  struct RapLess : public BoolParticleBaseFunctor {
    RapLess(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() < rapcut; }
    double rapcut;
  };

  /// Abs rapidity greater-than functor
  struct AbsRapGtr : public BoolParticleBaseFunctor {
    AbsRapGtr(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() > absrapcut; }
    double absrapcut;
  };

  /// Abs rapidity momentum less-than functor
  struct AbsRapLess : public BoolParticleBaseFunctor {
    AbsRapLess(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() < absrapcut; }
    double absrapcut;
  };


  /// Delta R (with respect to another 4-momentum, @a vec) greater-than functor
  struct DeltaRGtr : public BoolParticleBaseFunctor {
    DeltaRGtr(const FourMomentum& vec, double dr) : refvec(vec), drcut(dr) { }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec) > drcut; }
    FourMomentum refvec;
    double drcut;
  };

  /// Delta R (with respect to another 4-momentum, @a vec) less-than functor
  struct DeltaRLess : public BoolParticleBaseFunctor {
    DeltaRLess(const FourMomentum& vec, double dr) : refvec(vec), drcut(dr) { }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec) < drcut; }
    FourMomentum refvec;
    double drcut;
  };

  //@}


}

#endif
