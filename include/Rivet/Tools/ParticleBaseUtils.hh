#ifndef RIVET_PARTICLEBASEUTILS_HH
#define RIVET_PARTICLEBASEUTILS_HH

#include "Rivet/ParticleBase.hh"

namespace Rivet {



  /// @name ParticleBase classifier -> bool functors
  /// @todo Move to FourMomentum functions
  ///
  /// To be passed to any() or all() e.g. any(jets, DeltaRLess(electron, 0.4))
  //@{

  /// std::function instantiation for functors taking a ParticleBase and returning a bool
  using ParticleBaseSelector = function<bool(const ParticleBase&)>;
  /// std::function instantiation for functors taking two ParticleBase and returning a bool
  using ParticleBaseSorter = function<bool(const ParticleBase&, const ParticleBase&)>;


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
  using pTGtr = PtGtr;
  using ptGtr = PtGtr;

  /// Transverse momentum less-than functor
  struct PtLess : public BoolParticleBaseFunctor {
    PtLess(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() < ptcut; }
    double ptcut;
  };
  using pTLess = PtLess;
  using ptLess = PtLess;


  /// Pseudorapidity greater-than functor
  struct EtaGtr : public BoolParticleBaseFunctor {
    EtaGtr(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() > etacut; }
    double etacut;
  };
  using etaGtr = EtaGtr;

  /// Pseudorapidity momentum less-than functor
  struct EtaLess : public BoolParticleBaseFunctor {
    EtaLess(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() < etacut; }
    double etacut;
  };
  using etaLess = EtaLess;

  /// Abs pseudorapidity greater-than functor
  struct AbsEtaGtr : public BoolParticleBaseFunctor {
    AbsEtaGtr(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() > absetacut; }
    double absetacut;
  };
  using absEtaGtr = AbsEtaGtr;
  using absetaGtr = AbsEtaGtr;

  /// Abs pseudorapidity momentum less-than functor
  struct AbsEtaLess : public BoolParticleBaseFunctor {
    AbsEtaLess(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() < absetacut; }
    double absetacut;
  };
  using absEtaLess = AbsEtaLess;
  using absetaLess = AbsEtaLess;


  /// Rapidity greater-than functor
  struct RapGtr : public BoolParticleBaseFunctor {
    RapGtr(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() > rapcut; }
    double rapcut;
  };
  using rapGtr = RapGtr;

  /// Rapidity momentum less-than functor
  struct RapLess : public BoolParticleBaseFunctor {
    RapLess(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() < rapcut; }
    double rapcut;
  };
  using rapLess = RapLess;

  /// Abs rapidity greater-than functor
  struct AbsRapGtr : public BoolParticleBaseFunctor {
    AbsRapGtr(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() > absrapcut; }
    double absrapcut;
  };
  using absRapGtr = AbsRapGtr;
  using absrapGtr = AbsRapGtr;

  /// Abs rapidity momentum less-than functor
  struct AbsRapLess : public BoolParticleBaseFunctor {
    AbsRapLess(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() < absrapcut; }
    double absrapcut;
  };
  using absRapLess = AbsRapLess;
  using absrapLess = AbsRapLess;


  /// @f$ \Delta R @f$ (with respect to another 4-momentum, @a vec) greater-than functor
  struct DeltaRGtr : public BoolParticleBaseFunctor {
    DeltaRGtr(const ParticleBase& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec.mom()), drcut(dr), rapscheme(scheme) { }
    DeltaRGtr(const FourMomentum& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    DeltaRGtr(const Vector3& vec, double dr)
      : drcut(dr), rapscheme(PSEUDORAPIDITY) { refvec.setPx(vec.x()); refvec.setPy(vec.y()); refvec.setPz(vec.z()); }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec, rapscheme) > drcut; }
    FourMomentum refvec;
    double drcut;
    RapScheme rapscheme;
  };
  using deltaRGtr = DeltaRGtr;

  /// @f$ \Delta R @f$ (with respect to another 4-momentum, @a vec) less-than functor
  struct DeltaRLess : public BoolParticleBaseFunctor {
    DeltaRLess(const ParticleBase& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec.mom()), drcut(dr), rapscheme(scheme) { }
    DeltaRLess(const FourMomentum& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    DeltaRLess(const Vector3& vec, double dr)
      : drcut(dr), rapscheme(PSEUDORAPIDITY) { refvec.setPx(vec.x()); refvec.setPy(vec.y()); refvec.setPz(vec.z()); }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec, rapscheme) < drcut; }
    FourMomentum refvec;
    double drcut;
    RapScheme rapscheme;
  };
  using deltaRLess = DeltaRLess;


  /// @f$ |\Delta \phi| @f$ (with respect to another momentum, @a vec) greater-than functor
  struct DeltaPhiGtr : public BoolParticleBaseFunctor {
    DeltaPhiGtr(const ParticleBase& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiGtr(const FourMomentum& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiGtr(const Vector3& vec, double dphi)
      : refvec(vec), dphicut(dphi) { }
    bool operator()(const ParticleBase& p) const { return deltaPhi(p, refvec) > dphicut; }
    Vector3 refvec;
    double dphicut;
  };
  using deltaPhiGtr = DeltaPhiGtr;

  /// @f$ |\Delta \phi| @f$ (with respect to another momentum, @a vec) less-than functor
  struct DeltaPhiLess : public BoolParticleBaseFunctor {
    DeltaPhiLess(const ParticleBase& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiLess(const FourMomentum& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiLess(const Vector3& vec, double dphi)
      : refvec(vec), dphicut(dphi) { }
    bool operator()(const ParticleBase& p) const { return deltaPhi(p, refvec) < dphicut; }
    Vector3 refvec;
    double dphicut;
  };
  using deltaPhiLess = DeltaPhiLess;


  /// @f$ |\Delta \eta| @f$ (with respect to another momentum, @a vec) greater-than functor
  struct DeltaEtaGtr : public BoolParticleBaseFunctor {
    DeltaEtaGtr(const ParticleBase& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaGtr(const FourMomentum& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaGtr(const Vector3& vec, double deta)
      : refvec(vec), detacut(deta) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaEta(p, refvec)) > detacut; }
    Vector3 refvec;
    double detacut;
  };
  using deltaEtaGtr = DeltaEtaGtr;

  /// @f$ |\Delta \eta| @f$ (with respect to another momentum, @a vec) less-than functor
  struct DeltaEtaLess : public BoolParticleBaseFunctor {
    DeltaEtaLess(const ParticleBase& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaLess(const FourMomentum& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaLess(const Vector3& vec, double deta)
      : refvec(vec), detacut(deta) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaEta(p, refvec)) < detacut; }
    Vector3 refvec;
    double detacut;
  };
  using deltaEtaLess = DeltaEtaLess;


  /// @f$ |\Delta y| @f$ (with respect to another momentum, @a vec) greater-than functor
  struct DeltaRapGtr : public BoolParticleBaseFunctor {
    DeltaRapGtr(const ParticleBase& vec, double drap)
      : refvec(vec.mom()), drapcut(drap) { }
    DeltaRapGtr(const FourMomentum& vec, double drap)
      : refvec(vec), drapcut(drap) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaRap(p, refvec)) > drapcut; }
    FourMomentum refvec;
    double drapcut;
  };
  using deltaRapGtr = DeltaRapGtr;

  /// @f$ |\Delta y| @f$ (with respect to another momentum, @a vec) less-than functor
  struct DeltaRapLess : public BoolParticleBaseFunctor {
    DeltaRapLess(const ParticleBase& vec, double drap)
      : refvec(vec.mom()), drapcut(drap) { }
    DeltaRapLess(const FourMomentum& vec, double drap)
      : refvec(vec), drapcut(drap) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaRap(p, refvec)) < drapcut; }
    FourMomentum refvec;
    double drapcut;
  };
  using deltaRapLess = DeltaRapLess;

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

  /// Calculator of @f$ \Delta R @f$ with respect to a given momentum
  struct DeltaRWRT : public DoubleParticleBaseFunctor {
    DeltaRWRT(const ParticleBase& pb, RapScheme scheme=PSEUDORAPIDITY) : p(pb.mom()) {}
    DeltaRWRT(const FourMomentum& p4, RapScheme scheme=PSEUDORAPIDITY) : p(p4) {}
    DeltaRWRT(const Vector3& p3) : p(p3.mod(), p3.x(), p3.y(), p3.z()), rapscheme(PSEUDORAPIDITY) {}
    double operator()(const ParticleBase& pb) const { return deltaR(p, pb, rapscheme); }
    double operator()(const FourMomentum& p4) const { return deltaR(p, p4, rapscheme); }
    double operator()(const Vector3& p3) const { return deltaR(p, p3); }
    const FourMomentum p;
    RapScheme rapscheme;
  };
  using deltaRWRT = DeltaRWRT;

  /// Calculator of @f$ \Delta \phi @f$ with respect to a given momentum
  struct DeltaPhiWRT : public DoubleParticleBaseFunctor {
    DeltaPhiWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    DeltaPhiWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    DeltaPhiWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return deltaPhi(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaPhi(p, p4); }
    double operator()(const Vector3& p3) const { return deltaPhi(p, p3); }
    const Vector3 p;
  };
  using deltaPhiWRT = DeltaPhiWRT;

  /// Calculator of @f$ \Delta \eta @f$ with respect to a given momentum
  struct DeltaEtaWRT : public DoubleParticleBaseFunctor {
    DeltaEtaWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    DeltaEtaWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    DeltaEtaWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return deltaEta(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaEta(p, p4); }
    double operator()(const Vector3& p3) const { return deltaEta(p, p3); }
    const Vector3 p;
  };
  using deltaEtaWRT = DeltaEtaWRT;

  /// Calculator of @f$ |\Delta \eta| @f$ with respect to a given momentum
  struct AbsDeltaEtaWRT : public DoubleParticleBaseFunctor {
    AbsDeltaEtaWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    AbsDeltaEtaWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    AbsDeltaEtaWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return fabs(deltaEta(p, pb)); }
    double operator()(const FourMomentum& p4) const { return fabs(deltaEta(p, p4)); }
    double operator()(const Vector3& p3) const { return fabs(deltaEta(p, p3)); }
    const Vector3 p;
  };
  using absDeltaEtaWRT = AbsDeltaEtaWRT;

  /// Calculator of @f$ \Delta y @f$ with respect to a given momentum
  struct DeltaRapWRT : public DoubleParticleBaseFunctor {
    DeltaRapWRT(const ParticleBase& pb) : p(pb.mom()) {}
    DeltaRapWRT(const FourMomentum& p4) : p(p4) {}
    double operator()(const ParticleBase& pb) const { return deltaRap(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaRap(p, p4); }
    const FourMomentum p;
  };
  using deltaRapWRT = DeltaRapWRT;

  /// Calculator of @f$ |\Delta y| @f$ with respect to a given momentum
  struct AbsDeltaRapWRT : public DoubleParticleBaseFunctor {
    AbsDeltaRapWRT(const ParticleBase& pb) : p(pb.mom()) {}
    AbsDeltaRapWRT(const FourMomentum& p4) : p(p4) {}
    double operator()(const ParticleBase& pb) const { return fabs(deltaRap(p, pb)); }
    double operator()(const FourMomentum& p4) const { return fabs(deltaRap(p, p4)); }
    const FourMomentum p;
  };
  using absDeltaRapWRT = AbsDeltaRapWRT;

  //@}


  /// @name Non-PID particle properties, via unbound functions
  /// @todo Mostly move to functions on FourMomentum
  //@{

  /// Unbound function access to momentum
  inline FourMomentum mom(const ParticleBase& p) { return p.mom(); }

  /// Unbound function access to p3
  inline Vector3 p3(const ParticleBase& p) { return p.p3(); }

  /// Unbound function access to pTvec
  inline Vector3 pTvec(const ParticleBase& p) { return p.pTvec(); }

  /// Unbound function access to p
  inline double p(const ParticleBase& p) { return p.p(); }

  /// Unbound function access to pT
  inline double pT(const ParticleBase& p) { return p.pT(); }

  /// Unbound function access to ET
  inline double Et(const ParticleBase& p) { return p.Et(); }

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
