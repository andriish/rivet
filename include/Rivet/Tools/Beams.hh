// -*- C++ -*-
#ifndef RIVET_Beams_HH
#define RIVET_Beams_HH

#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {


  /// @defgroup beam_functions Standalone beam kinematics functions
  /// @{

  /// Get beam particles from an event
  ParticlePair beams(const Event& e);


  /// Check if the given beam pair is valid
  bool validBeams(const ParticlePair& beams);

  /// Check if the event's beam pair is valid
  inline bool validBeams(const Event& e) {
    return validBeams(beams(e));
  }


  /// Get beam particle IDs from a pair of Particles
  inline PdgIdPair beamIDs(const ParticlePair& beams) { return pids(beams); }

  /// Get beam particle IDs from an event
  inline PdgIdPair beamIDs(const Event& e) { return beamIDs(beams(e)); }


  /// Get beam particle energies from a pair of Particles
  inline PdgIdPair beamEnergies(const ParticlePair& beams) { return energies(beams); }

  /// Get beam particle energies from an event
  inline PdgIdPair beamEnergies(const Event& e) { return beamEnergies(beams(e)); }


  /// @brief Get beam centre-of-mass energy from a pair of beam energies
  ///
  /// For want of more complete information, this function assumes massless
  /// particles colliding head-on.
  double sqrtS(double ea, double eb);

  /// @brief Get beam centre-of-mass energy from a pair of beam energies
  ///
  /// For want of more complete information, this function assumes massless
  /// particles colliding head-on.
  inline double sqrtS(const pair<double,double>& energies) { return sqrtS(energies.first, energies.second); }

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
  inline FourMomentum cmsBoostVec(const FourMomentum& pa, const FourMomentum& pb) {
    return pa + pb;
  }

  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of Particles
  inline FourMomentum cmsBoostVec(const ParticlePair& beams) {
    return cmsBoostVec(beams.first, beams.second);
  }

  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of beam momenta
  FourMomentum acmsBoostVec(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of Particles
  FourMomentum acmsBoostVec(const ParticlePair& beams);


  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of beam momenta
  Vector3 cmsBetaVec(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of Particles
  inline Vector3 cmsBetaVec(const ParticlePair& beams) {
    return cmsBetaVec(beams.first, beams.second);
  }


  /// Get the Lorentz boost to the per-nucleon beam centre-of-mass system (ACMS) from a pair of beam momenta
  /// @note Uses a nominal nucleon mass of 0.939 GeV to convert masses to A
  Vector3 acmsBetaVec(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz boost to the per-nucleon beam centre-of-mass system (ACMS) from a pair of Particles
  /// @note Uses the sum of nuclear mass numbers A for each beam
  Vector3 acmsBetaVec(const ParticlePair& beams);


  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of beam momenta
  Vector3 cmsGammaVec(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz boost to the beam centre-of-mass system (CMS) from a pair of Particles
  inline Vector3 cmsGammaVec(const ParticlePair& beams) {
    return cmsGammaVec(beams.first, beams.second);
  }


  /// Get the Lorentz boost to the per-nucleon beam centre-of-mass system (ACMS) from a pair of beam momenta
  /// @note Uses a nominal nucleon mass of 0.939 GeV to convert masses to A
  Vector3 acmsGammaVec(const FourMomentum& pa, const FourMomentum& pb);

  /// Get the Lorentz boost to the per-nucleon beam centre-of-mass system (ACMS) from a pair of Particles
  /// @note Uses the sum of nuclear mass numbers A for each beam
  Vector3 acmsGammaVec(const ParticlePair& beams);


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

  /// @}


}

#endif
