// -*- C++ -*-
#ifndef RIVET_Beam_HH
#define RIVET_Beam_HH

#include "Rivet/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {


  /// @name Standalone beam kinematics functions
  //@{

  /// Get beam particles from an event
  ParticlePair beams(const Event& e);

  /// Get beam particle IDs from a pair of Particles
  /// @deprecated Use pids(beams)
  inline PdgIdPair beamIds(const ParticlePair& beams) { return pids(beams); }

  /// Get beam particle IDs from an event
  /// @deprecated Use pids(e.beams())
  inline PdgIdPair beamIds(const Event& e) { return pids(beams(e)); }


  /// Get beam centre-of-mass energy from a pair of beam momenta
  double sqrtS(const FourMomentum& pa, const FourMomentum& pb);

  /// Get beam centre-of-mass energy from a pair of Particles
  inline double sqrtS(const ParticlePair& beams) {
    return sqrtS(beams.first.momentum(), beams.second.momentum());
  }

  /// Get beam centre-of-mass energy from an Event
  inline double sqrtS(const Event& e) { return sqrtS(beams(e)); }


  /// Get per-nucleon beam centre-of-mass energy from a pair of beam momenta
  /// @note Uses a nominal nucleon mass of 0.939 GeV to convert masses to A
  double asqrtS(const FourMomentum& pa, const FourMomentum& pb);

  /// Get per-nucleon beam centre-of-mass energy from a pair of Particles
  /// @note Uses the sum of nuclear mass numbers A for each beam
  double asqrtS(const ParticlePair& beams);

  /// Get per-nucleon beam centre-of-mass energy from an Event
  /// @note Uses the sum of nuclear mass numbers A for each beam
  inline double asqrtS(const Event& e) { return asqrtS(beams(e)); }


  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of beam momenta
  Vector3 cmsBoostBetaVec(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of Particles
  inline Vector3 cmsBoostBetaVec(const ParticlePair& beams) {
    return cmsBoostBetaVec(beams.first, beams.second);
  }


  /// Get the Lorentz boost to the per-nucleon beam centre-of-mass system (ACMS) from a pair of beam momenta
  /// @note Uses a nominal nucleon mass of 0.939 GeV to convert masses to A
  Vector3 acmsBoostBetaVec(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz boost to the per-nucleon beam centre-of-mass system (ACMS) from a pair of Particles
  /// @note Uses the sum of nuclear mass numbers A for each beam
  Vector3 acmsBoostBetaVec(const ParticlePair& beams);


  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of beam momenta
  Vector3 cmsBoostGammaVec(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of Particles
  inline Vector3 cmsBoostGammaVec(const ParticlePair& beams) {
    return cmsBoostGammaVec(beams.first, beams.second);
  }


  /// Get the Lorentz boost to the per-nucleon beam centre-of-mass system (ACMS) from a pair of beam momenta
  /// @note Uses a nominal nucleon mass of 0.939 GeV to convert masses to A
  Vector3 acmsBoostGammaVec(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz boost to the per-nucleon beam centre-of-mass system (ACMS) from a pair of Particles
  /// @note Uses the sum of nuclear mass numbers A for each beam
  Vector3 acmsBoostGammaVec(const ParticlePair& beams);


  /// Get the Lorentz transformation to the beam centre-of-mass system (CMS) from a pair of beam momenta
  LorentzTransform cmsTransform(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz transformation to the beam centre-of-mass system (CMS) from a pair of Particles
  inline LorentzTransform cmsTransform(const ParticlePair& beams) {
    return cmsTransform(beams.first, beams.second);
  }


  /// Get the Lorentz transformation to the per-nucleon beam centre-of-mass system (CMS) from a pair of beam momenta
  /// @note Uses a nominal nucleon mass of 0.939 GeV to convert masses to A
  LorentzTransform acmsTransform(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz transformation to the per-nucleon beam centre-of-mass system (CMS) from a pair of Particles
  /// @note Uses the sum of nuclear mass numbers A for each beam
  LorentzTransform acmsTransform(const ParticlePair& beams);

  //@}




  /// @brief Project out the incoming beams
  class Beam : public Projection {
  public:

    /// Default (and only) constructor
    Beam() { setName("Beam"); }

    /// Clone on the heap
    DEFAULT_RIVET_PROJ_CLONE(Beam);


    /// @name Beam particles and kinematics
    //@{

    /// The pair of beam particles in the current collision
    const ParticlePair& beams() const { return _theBeams; }

    /// The pair of beam particle PDG codes in the current collision
    /// @deprecated Use pids(beams())
    PdgIdPair beamIds() const { return pids(beams()); }

    /// Get centre of mass energy, \f$ \sqrt{s} \f$
    double sqrtS() const { return Rivet::sqrtS(beams()); }

    /// Get the Lorentz boost to the beam centre-of-mass
    Vector3 cmsBoostBeta() const { return Rivet::cmsBoostBetaVec(beams()); }

    /// Get the Lorentz boost to the beam centre-of-mass
    Vector3 cmsBoostGamma() const { return Rivet::cmsBoostGammaVec(beams()); }

    /// Get the Lorentz transform to the beam centre-of-mass
    LorentzTransform cmsTransform() const { return Rivet::cmsTransform(beams()); }

    //@}


    /// @name Per-nucleon beam kinematics
    //@{

    /// Get per-nucleon centre of mass energy, \f$ \sqrt{s}/(A_1 + A_2) \f$
    double asqrtS() const { return Rivet::asqrtS(beams()); }

    /// Get the Lorentz boost to the per-nucleon beam centre-of-mass
    Vector3 acmsBoostBeta() const { return Rivet::acmsBoostBetaVec(beams()); }

    /// Get the Lorentz boost to the per-nucleon beam centre-of-mass
    Vector3 acmsBoostGamma() const { return Rivet::acmsBoostGammaVec(beams()); }

    /// Get the Lorentz transform to the per-nucleon beam centre-of-mass
    LorentzTransform acmsTransform() const { return Rivet::acmsTransform(beams()); }

    //@}


    /// Get the beam interaction primary vertex (PV) position
    FourVector pv() const;


  public:

    /// Project on to the Event
    virtual void project(const Event& e);


  private:

    /// Compare with other projections -- it's always the same, since there are no params
    virtual int compare(const Projection& UNUSED(p)) const { return EQUIVALENT; }

    /// The beam particles in the current collision
    ParticlePair _theBeams;

  };


}

#endif
