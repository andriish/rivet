// -*- C++ -*-
#ifndef RIVET_Beam_HH
#define RIVET_Beam_HH

#include "Rivet/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  /// @name Standalone beam kinematics functions
  //@{

  /// Get beam particles from an event
  ParticlePair beams(const Event& e);


  /// Get beam particle IDs from a pair of Particles
  inline PdgIdPair beamIds(const ParticlePair& beams) {
    return make_pair(beams.first.pid(), beams.second.pid());
  }

  /// Get beam particle IDs from an event
  inline PdgIdPair beamIds(const Event& e) {
    return beamIds(beams(e));
  }


  /// Get beam centre-of-mass energy from a pair of beam momenta
  double sqrtS(const FourMomentum& pa, const FourMomentum& pb);

  /// Get beam centre-of-mass energy from a pair of Particles
  inline double sqrtS(const ParticlePair& beams) {
    return sqrtS(beams.first.momentum(), beams.second.momentum());
  }

  /// Get beam centre-of-mass energy from an Event
  inline double sqrtS(const Event& e) {
    return sqrtS(beams(e));
  }


  /// Get per-nucleon beam centre-of-mass energy from a pair of beam momenta
  /// @note Uses a nominal nucleon mass of 0.939 GeV to convert masses to A
  double asqrtS(const FourMomentum& pa, const FourMomentum& pb);

  /// Get per-nucleon beam centre-of-mass energy from a pair of Particles
  /// @note Uses the sum of nuclear mass numbers A for each beam
  double asqrtS(const ParticlePair& beams);

  /// Get per-nucleon beam centre-of-mass energy from an Event
  inline double asqrtS(const Event& e) {
    return asqrtS(beams(e));
  }


  /// Get an active Lorentz boost to the beam centre-of-mass from a pair of beam momenta
  Vector3 comBoost(const FourMomentum& pa, const FourMomentum& pb);

  /// Get an active Lorentz boost to the beam centre-of-mass from a pair of Particles
  inline Vector3 comBoost(const ParticlePair& beams) {
    return comBoost(beams.first, beams.second);
  }

  /// Get an active Lorentz boost to the beam centre-of-mass from an Event
  inline Vector3 comBoost(const Event& e) {
    return comBoost(beams(e));
  }

  //@}




  /// @brief Project out the incoming beams
  class Beam : public Projection {
  public:

    /// The default constructor.
    Beam() {
      setName("Beam");
    }

    /// Clone on the heap.
    virtual unique_ptr<Projection> clone() const {
      return unique_ptr<Projection>(new Beam(*this));
    }


  public:

    /// The pair of beam particles in the current collision.
    const ParticlePair& beams() const {
      return _theBeams;
    }

    /// The pair of beam particle PDG codes in the current collision.
    const PdgIdPair beamIds() const {
      return Rivet::beamIds(beams());
    }

    /// Get centre of mass energy, \f$ \sqrt{s} \f$.
    double sqrtS() const;

    /// Get the beam interaction primary vertex (PV) position.
    FourVector pv() const;


  public:

    /// Project on to the Event
    virtual void project(const Event& e);


  protected:

    /// Compare with other projections.
    virtual int compare(const Projection& UNUSED(p)) const {
      return EQUIVALENT;
    }


  private:

    /// The beam particles in the current collision
    ParticlePair _theBeams;

  };


}

#endif
