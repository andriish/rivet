// -*- C++ -*-
#ifndef RIVET_Jet_HH
#define RIVET_Jet_HH

#include "Rivet/Rivet.hh"
#include <numeric>

namespace Rivet {


  /// @brief Representation of a clustered jet of particles.
  class Jet : public ParticleBase {
  public:

    /// @name Constructors
    //@{

    Jet() : ParticleBase() { clear(); }

    /// Set all the jet data, with full particle information.
    Jet(const vector<Particle>& particles, const FourMomentum& pjet)
      : ParticleBase() {
      setState(particles, pjet);
    }

    /// @todo Add a constructor from PseudoJet

    //@}


    /// @name Access jet constituents
    //@{

    /// Number of particles in this jet.
    size_t size() const { return _particles.size(); }

    /// Get the particles in this jet.
    vector<Particle>& particles() { return _particles; }

    /// Get the particles in this jet (const version)
    const vector<Particle>& particles() const { return _particles; }

    /// Check whether this jet contains a particular particle.
    bool containsParticle(const Particle& particle) const;

    /// Check whether this jet contains a certain particle type.
    bool containsParticleId(PdgId pid) const;

    /// Check whether this jet contains at least one of certain particle types.
    bool containsParticleId(const vector<PdgId>& pids) const;

    /// Check whether this jet contains a charm-flavoured hadron (or decay products from one).
    bool containsCharm() const;

    /// Check whether this jet contains a bottom-flavoured hadron (or decay products from one).
    bool containsBottom() const;

    //@}


    /// @name Access additional effective jet 4-vector properties
    //@{

    /// Get equivalent single momentum four-vector.
    const FourMomentum& momentum() const { return _momentum; }

    /// Get the total energy of this jet.
    double totalEnergy() const { return momentum().E(); }

    /// Get the energy carried in this jet by neutral particles.
    double neutralEnergy() const;

    /// Get the energy carried in this jet by hadrons.
    double hadronicEnergy() const;

    //@}


    /// @name Set the jet constituents and properties
    //@{

    /// Set all the jet data, with full particle information.
    Jet& setState(const vector<Particle>& particles, const FourMomentum& pjet);

    /// Set the effective 4-momentum of the jet.
    Jet& setMomentum(const FourMomentum& momentum);

    /// Set the particles collection with full particle information.
    Jet& setParticles(const vector<Particle>& particles);

    /// Reset this jet as empty.
    Jet& clear();

    //@}


  private:

    /// Full particle information including tracks, ID etc
    Particles _particles;

    /// Effective jet 4-vector
    FourMomentum _momentum;

    /// @todo Add a FJ3 PseudoJet member to unify PseudoJet and Jet

  };


  /// Typedef for a collection of Jet objects.
  typedef vector<Jet> Jets;


  /// @name deltaR, deltaEta, deltaPhi functions specifically for Jet arguments
  //@{

  inline double deltaR(const Jet& j1, const Jet& j2,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(j1.momentum(), j2.momentum(), scheme);
  }

  inline double deltaR(const Jet& j, const Particle& p,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(j.momentum(), p.momentum(), scheme);
  }

  inline double deltaR(const Particle& p, const Jet& j,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(p.momentum(), j.momentum(), scheme);
  }

  inline double deltaR(const Jet& j, const FourMomentum& v,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(j.momentum(), v, scheme);
  }

  inline double deltaR(const Jet& j, const FourVector& v,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(j.momentum(), v, scheme);
  }

  inline double deltaR(const Jet& j, const Vector3& v) {
    return deltaR(j.momentum(), v);
  }

  inline double deltaR(const Jet& j, double eta, double phi) {
    return deltaR(j.momentum(), eta, phi);
  }

  inline double deltaR(const FourMomentum& v, const Jet& j,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(v, j.momentum(), scheme);
  }

  inline double deltaR(const FourVector& v, const Jet& j,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(v, j.momentum(), scheme);
  }

  inline double deltaR(const Vector3& v, const Jet& j) {
    return deltaR(v, j.momentum());
  }

  inline double deltaR(double eta, double phi, const Jet& j) {
    return deltaR(eta, phi, j.momentum());
  }


  inline double deltaPhi(const Jet& j1, const Jet& j2) {
    return deltaPhi(j1.momentum(), j2.momentum());
  }

  inline double deltaPhi(const Jet& j, const Particle& p) {
    return deltaPhi(j.momentum(), p.momentum());
  }

  inline double deltaPhi(const Particle& p, const Jet& j) {
    return deltaPhi(p.momentum(), j.momentum());
  }

  inline double deltaPhi(const Jet& j, const FourMomentum& v) {
    return deltaPhi(j.momentum(), v);
  }

  inline double deltaPhi(const Jet& j, const FourVector& v) {
    return deltaPhi(j.momentum(), v);
  }

  inline double deltaPhi(const Jet& j, const Vector3& v) {
    return deltaPhi(j.momentum(), v);
  }

  inline double deltaPhi(const Jet& j, double phi) {
    return deltaPhi(j.momentum(), phi);
  }

  inline double deltaPhi(const FourMomentum& v, const Jet& j) {
    return deltaPhi(v, j.momentum());
  }

  inline double deltaPhi(const FourVector& v, const Jet& j) {
    return deltaPhi(v, j.momentum());
  }

  inline double deltaPhi(const Vector3& v, const Jet& j) {
    return deltaPhi(v, j.momentum());
  }

  inline double deltaPhi(double phi, const Jet& j) {
    return deltaPhi(phi, j.momentum());
  }


  inline double deltaEta(const Jet& j1, const Jet& j2) {
    return deltaEta(j1.momentum(), j2.momentum());
  }

  inline double deltaEta(const Jet& j, const Particle& p) {
    return deltaEta(j.momentum(), p.momentum());
  }

  inline double deltaEta(const Particle& p, const Jet& j) {
    return deltaEta(p.momentum(), j.momentum());
  }

  inline double deltaEta(const Jet& j, const FourMomentum& v) {
    return deltaEta(j.momentum(), v);
  }

  inline double deltaEta(const Jet& j, const FourVector& v) {
    return deltaEta(j.momentum(), v);
  }

  inline double deltaEta(const Jet& j, const Vector3& v) {
    return deltaEta(j.momentum(), v);
  }

  inline double deltaEta(const Jet& j, double eta) {
    return deltaEta(j.momentum(), eta);
  }

  inline double deltaEta(const FourMomentum& v, const Jet& j) {
    return deltaEta(v, j.momentum());
  }

  inline double deltaEta(const FourVector& v, const Jet& j) {
    return deltaEta(v, j.momentum());
  }

  inline double deltaEta(const Vector3& v, const Jet& j) {
    return deltaEta(v, j.momentum());
  }

  inline double deltaEta(double eta, const Jet& j) {
    return deltaEta(eta, j.momentum());
  }

  //@}


}

#endif
