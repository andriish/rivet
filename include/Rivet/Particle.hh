// -*- C++ -*-
#ifndef RIVET_Particle_HH
#define RIVET_Particle_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/ParticleBase.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  /// Representation of particles from a HepMC::GenEvent.
  class Particle : public ParticleBase {
  public:

    /// Default constructor.
    /// @deprecated A particle without info is useless. This only exists to keep STL containers happy.
    Particle()
      : ParticleBase(),
        _original(0), _id(0), _momentum()
    { }

    /// Constructor without GenParticle.
    Particle(PdgId pid, const FourMomentum& mom)
      : ParticleBase(),
        _original(0), _id(pid), _momentum(mom)
    { }

    /// Constructor from a HepMC GenParticle.
    Particle(const GenParticle& gp)
      : ParticleBase(),
        _original(&gp), _id(gp.pdg_id()),
        _momentum(gp.momentum())
    { }

    /// Constructor from a HepMC GenParticle pointer.
    Particle(const GenParticle* gp)
      : ParticleBase(),
        _original(gp), _id(gp->pdg_id()),
        _momentum(gp->momentum())
    { }


  public:

    /// @name Basic particle specific properties
    //@{

    /// Get a const reference to the original GenParticle.
    const GenParticle* genParticle() const {
      return _original;
    }

    /// Cast operator for conversion to GenParticle*
    operator const GenParticle* () const { return genParticle(); }

    /// The momentum.
    const FourMomentum& momentum() const {
      return _momentum;
    }
    /// Set the momentum.
    Particle& setMomentum(const FourMomentum& momentum) {
      _momentum = momentum;
      return *this;
    }

    /// This Particle's PDG ID code.
    PdgId pdgId() const {
      return _id;
    }

    //@}

    /// @name Particle species properties
    //@{

    /// The charge of this Particle.
    double charge() const {
      return PID::charge(pdgId());
    }
    /// Three times the charge of this Particle (i.e. integer multiple of smallest quark charge).
    int threeCharge() const {
      return PID::threeCharge(pdgId());
    }

    /// Is this a hadron?
    bool isHadron() const { return PID::isHadron(pdgId()); }

    /// Is this a meson?
    bool isMeson() const { return PID::isMeson(pdgId()); }

    /// Is this a baryon?
    bool isBaryon() const { return PID::isBaryon(pdgId()); }

    /// Is this a lepton?
    bool isLepton() const { return PID::isLepton(pdgId()); }

    /// Is this a neutrino?
    bool isNeutrino() const { return PID::isNeutrino(pdgId()); }

    /// Does this (hadron) contain a b quark?
    bool hasBottom() const { return PID::hasBottom(pdgId()); }

    /// Does this (hadron) contain a c quark?
    bool hasCharm() const { return PID::hasCharm(pdgId()); }

    // /// Does this (hadron) contain an s quark?
    // bool hasStrange() const { return PID::hasStrange(pdgId()); }

    //@}


    /// @name Ancestry properties
    //@{

    /// Check whether a given PID is found in the GenParticle's ancestor list
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasAncestor(PdgId pdg_id) const;

    /// @brief Determine whether the particle is from a hadron or tau decay
    ///
    /// Specifically, walk up the ancestor chain until a status 2 hadron or
    /// tau is found, if at all.
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromDecay() const;

    /// @brief Determine whether the particle is from a b-hadron decay
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromBottom() const;

    /// @brief Determine whether the particle is from a c-hadron decay
    ///
    /// @note If a hadron contains b and c quarks it is considered a bottom
    /// hadron and NOT a charm hadron.
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromCharm() const;

    // /// @brief Determine whether the particle is from a s-hadron decay
    // ///
    // /// @note If a hadron contains b or c quarks as well as strange it is
    // /// considered a b or c hadron, but NOT a strange hadron.
    // ///
    // /// @note This question is valid in MC, but may not be perfectly answerable
    // /// experimentally -- use this function with care when replicating
    // /// experimental analyses!
    // bool fromStrange() const;

    /// @brief Determine whether the particle is from a tau decay
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromTau() const;

    //@}


  private:

    /// A pointer to the original GenParticle from which this Particle is projected.
    const GenParticle* _original;

    /// The PDG ID code for this Particle.
    PdgId _id;

    /// The momentum of this projection of the Particle.
    FourMomentum _momentum;

    /// @todo Also store production and decay positions and make them available.

  };


  /// @name String representation
  //@{

  /// Print a ParticlePair as a string.
  inline std::string toString(const ParticlePair& pair) {
    stringstream out;
    out << "["
        << PID::toParticleName(pair.first.pdgId()) << " @ "
        << pair.first.momentum().E()/GeV << " GeV, "
        << PID::toParticleName(pair.second.pdgId()) << " @ "
        << pair.second.momentum().E()/GeV << " GeV]";
    return out.str();
  }

  /// Allow ParticlePair to be passed to an ostream.
  inline std::ostream& operator<<(std::ostream& os, const ParticlePair& pp) {
    os << toString(pp);
    return os;
  }

  //@}


  /// @name deltaR, deltaEta, deltaPhi functions specifically for Particle arguments
  //@{

  inline double deltaR(const Particle& p1, const Particle& p2,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(p1.momentum(), p2.momentum(), scheme);
  }

  inline double deltaR(const Particle& p, const FourMomentum& v,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(p.momentum(), v, scheme);
  }

  inline double deltaR(const Particle& p, const FourVector& v,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(p.momentum(), v, scheme);
  }

  inline double deltaR(const Particle& p, const Vector3& v) {
    return deltaR(p.momentum(), v);
  }

  inline double deltaR(const Particle& p, double eta, double phi) {
    return deltaR(p.momentum(), eta, phi);
  }

  inline double deltaR(const FourMomentum& v, const Particle& p,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(v, p.momentum(), scheme);
  }

  inline double deltaR(const FourVector& v, const Particle& p,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(v, p.momentum(), scheme);
  }

  inline double deltaR(const Vector3& v, const Particle& p) {
    return deltaR(v, p.momentum());
  }

  inline double deltaR(double eta, double phi, const Particle& p) {
    return deltaR(eta, phi, p.momentum());
  }


  inline double deltaPhi(const Particle& p1, const Particle& p2) {
    return deltaPhi(p1.momentum(), p2.momentum());
  }

  inline double deltaPhi(const Particle& p, const FourMomentum& v) {
    return deltaPhi(p.momentum(), v);
  }

  inline double deltaPhi(const Particle& p, const FourVector& v) {
    return deltaPhi(p.momentum(), v);
  }

  inline double deltaPhi(const Particle& p, const Vector3& v) {
    return deltaPhi(p.momentum(), v);
  }

  inline double deltaPhi(const Particle& p, double phi) {
    return deltaPhi(p.momentum(), phi);
  }

  inline double deltaPhi(const FourMomentum& v, const Particle& p) {
    return deltaPhi(v, p.momentum());
  }

  inline double deltaPhi(const FourVector& v, const Particle& p) {
    return deltaPhi(v, p.momentum());
  }

  inline double deltaPhi(const Vector3& v, const Particle& p) {
    return deltaPhi(v, p.momentum());
  }

  inline double deltaPhi(double phi, const Particle& p) {
    return deltaPhi(phi, p.momentum());
  }


  inline double deltaEta(const Particle& p1, const Particle& p2) {
    return deltaEta(p1.momentum(), p2.momentum());
  }

  inline double deltaEta(const Particle& p, const FourMomentum& v) {
    return deltaEta(p.momentum(), v);
  }

  inline double deltaEta(const Particle& p, const FourVector& v) {
    return deltaEta(p.momentum(), v);
  }

  inline double deltaEta(const Particle& p, const Vector3& v) {
    return deltaEta(p.momentum(), v);
  }

  inline double deltaEta(const Particle& p, double eta) {
    return deltaEta(p.momentum(), eta);
  }

  inline double deltaEta(const FourMomentum& v, const Particle& p) {
    return deltaEta(v, p.momentum());
  }

  inline double deltaEta(const FourVector& v, const Particle& p) {
    return deltaEta(v, p.momentum());
  }

  inline double deltaEta(const Vector3& v, const Particle& p) {
    return deltaEta(v, p.momentum());
  }

  inline double deltaEta(double eta, const Particle& p) {
    return deltaEta(eta, p.momentum());
  }

  //@}


}

#endif
