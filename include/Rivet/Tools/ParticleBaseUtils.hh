#ifndef RIVET_PARTICLEBASEUTILS_HH
#define RIVET_PARTICLEBASEUTILS_HH

#include "Rivet/ParticleBase.hh"

namespace Rivet {


  /// @name ParticleBase classifier -> bool functors
  /// @todo Move to FourMomentum functions
  ///
  /// To be passed to any() or all() e.g. any(jets, DeltaRLess(electron, 0.4))
  //@{

  /// Base type for Particle -> bool functors
  struct BoolParticleBaseFunctor {
    virtual bool operator()(const ParticleBase& p) const = 0;
  };

  /// Transverse momentum greater-than functor
  struct pTGtr : public BoolParticleBaseFunctor {
    pTGtr(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() > ptcut; }
    double ptcut;
  };
  using PtGtr = pTGtr;

  /// Transverse momentum less-than functor
  struct pTLess : public BoolParticleBaseFunctor {
    pTLess(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() < ptcut; }
    double ptcut;
  };
  using PtLess = pTLess;


  /// Pseudorapidity greater-than functor
  struct etaGtr : public BoolParticleBaseFunctor {
    etaGtr(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() > etacut; }
    double etacut;
  };
  using EtaGtr = etaGtr;

  /// Pseudorapidity momentum less-than functor
  struct etaLess : public BoolParticleBaseFunctor {
    etaLess(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() < etacut; }
    double etacut;
  };
  using EtaLess = etaLess;

  /// Abs pseudorapidity greater-than functor
  struct absetaGtr : public BoolParticleBaseFunctor {
    absetaGtr(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() > absetacut; }
    double absetacut;
  };
  using AbsEtaGtr = etaGtr;

  /// Abs pseudorapidity momentum less-than functor
  struct absetaLess : public BoolParticleBaseFunctor {
    absetaLess(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() < absetacut; }
    double absetacut;
  };
  using AbsEtaLess = absetaLess;


  /// Rapidity greater-than functor
  struct rapGtr : public BoolParticleBaseFunctor {
    rapGtr(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() > rapcut; }
    double rapcut;
  };
  using RapGtr = rapGtr;

  /// Rapidity momentum less-than functor
  struct rapLess : public BoolParticleBaseFunctor {
    rapLess(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() < rapcut; }
    double rapcut;
  };
  using RapLess = rapLess;

  /// Abs rapidity greater-than functor
  struct absrapGtr : public BoolParticleBaseFunctor {
    absrapGtr(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() > absrapcut; }
    double absrapcut;
  };
  using AbsRapGtr = absrapGtr;

  /// Abs rapidity momentum less-than functor
  struct absrapLess : public BoolParticleBaseFunctor {
    absrapLess(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() < absrapcut; }
    double absrapcut;
  };
  using AbsRapLess = absrapLess;


  /// Delta R (with respect to another 4-momentum, @a vec) greater-than functor
  struct deltaRGtr : public BoolParticleBaseFunctor {
    deltaRGtr(const FourMomentum& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec, rapscheme) > drcut; }
    FourMomentum refvec;
    double drcut;
    RapScheme rapscheme;
  };
  using DeltaRGtr = deltaRGtr;

  /// Delta R (with respect to another 4-momentum, @a vec) less-than functor
  struct deltaRLess : public BoolParticleBaseFunctor {
    deltaRLess(const FourMomentum& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec, rapscheme) < drcut; }
    FourMomentum refvec;
    double drcut;
    RapScheme rapscheme;
  };
  using DeltaRLess = deltaRLess;

  //@}


  /// @name ParticleBase comparison -> double functors
  /// @todo Move to FourMomentum functions
  ///
  /// To be passed to transform()any(jets, DeltaRLess(electron, 0.4))
  //@{

  /// Base type for Particle -> double functors
  struct DoubleParticleBaseFunctor {
    virtual double operator()(const ParticleBase& p) const = 0;
  };

  struct deltaRWRT : public DoubleParticleBaseFunctor {
    deltaRWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    deltaRWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    deltaRWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return deltaR(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaR(p, p4); }
    double operator()(const Vector3& p3) const { return deltaR(p, p3); }
    const Vector3 p;
    RapScheme rapscheme;
  };
  using DeltaRWRT = deltaRWRT;

  struct deltaPhiWRT : public DoubleParticleBaseFunctor {
    deltaPhiWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    deltaPhiWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    deltaPhiWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return deltaPhi(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaPhi(p, p4); }
    double operator()(const Vector3& p3) const { return deltaPhi(p, p3); }
    const Vector3 p;
  };
  using DeltaPhiWRT = deltaPhiWRT;

  struct deltaEtaWRT : public DoubleParticleBaseFunctor {
    deltaEtaWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    deltaEtaWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    deltaEtaWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return deltaEta(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaEta(p, p4); }
    double operator()(const Vector3& p3) const { return deltaEta(p, p3); }
    const Vector3 p;
  };
  using DeltaEtaWRT = deltaEtaWRT;

  struct absDeltaEtaWRT : public DoubleParticleBaseFunctor {
    absDeltaEtaWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    absDeltaEtaWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    absDeltaEtaWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return fabs(deltaEta(p, pb)); }
    double operator()(const FourMomentum& p4) const { return fabs(deltaEta(p, p4)); }
    double operator()(const Vector3& p3) const { return fabs(deltaEta(p, p3)); }
    const Vector3 p;
  };
  using absDeltaEtaWRT = absDeltaEtaWRT;

  struct deltaRapWRT : public DoubleParticleBaseFunctor {
    deltaRapWRT(const ParticleBase& pb) : p(pb.mom()) {}
    deltaRapWRT(const FourMomentum& p4) : p(p4) {}
    double operator()(const ParticleBase& pb) const { return deltaRap(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaRap(p, p4); }
    const FourMomentum p;
  };
  using DeltaRapWRT = deltaRapWRT;

  struct absDeltaRapWRT : public DoubleParticleBaseFunctor {
    absDeltaRapWRT(const ParticleBase& pb) : p(pb.mom()) {}
    absDeltaRapWRT(const FourMomentum& p4) : p(p4) {}
    double operator()(const ParticleBase& pb) const { return fabs(deltaRap(p, pb)); }
    double operator()(const FourMomentum& p4) const { return fabs(deltaRap(p, p4)); }
    const FourMomentum p;
  };
  using absDeltaRapWRT = absDeltaRapWRT;

  //@}


  /// @name Non-PID particle properties, via unbound functions
  /// @todo Move to FourMomentum functions
  //@{

  /// Unbound function access to momentum
  inline FourMomentum mom(const ParticleBase& p) { return p.mom(); }

  /// Unbound function access to p3
  inline Vector3 p3(const ParticleBase& p) { return p.p3(); }

  /// Unbound function access to p
  inline double p(const ParticleBase& p) { return p.p(); }

  /// Unbound function access to pT
  inline double pT(const ParticleBase& p) { return p.pT(); }

  /// Unbound function access to eta
  inline double eta(const ParticleBase& p) { return p.eta(); }

  /// Unbound function access to abseta
  inline double abseta(const ParticleBase& p) { return p.abseta(); }

  /// Unbound function access to rapidity
  inline double rap(const ParticleBase& p) { return p.rap(); }

  /// Unbound function access to abs rapidity
  inline double absrap(const ParticleBase& p) { return p.absrap(); }

  //@}


}

#endif
