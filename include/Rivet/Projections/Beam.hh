// -*- C++ -*-
#ifndef RIVET_Beam_HH
#define RIVET_Beam_HH

#include "Rivet/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Tools/Beams.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {


  /// @brief Project out the incoming beams
  class Beam : public Projection {
  public:

    /// Default (and only) constructor
    Beam() { setName("Beam"); }

    /// Clone on the heap
    DEFAULT_RIVET_PROJ_CLONE(Beam);


    /// @name Beam particles and kinematics
    /// @{

    /// The pair of beam particles in the current collision
    const ParticlePair& beams() const { return _theBeams; }

    /// The pair of beam particle PDG codes in the current collision
    /// @deprecated Use pids(beams())
    PdgIdPair beamIDs() const { return pids(beams()); }

    /// Get centre of mass energy, \f$ \sqrt{s} \f$
    double sqrtS() const { return Rivet::sqrtS(beams()); }

    /// Get the Lorentz boost to the beam centre-of-mass
    FourMomentum cmsBoostVec() const { return Rivet::cmsBoostVec(beams()); }

    /// Get the Lorentz transform to the beam centre-of-mass
    LorentzTransform cmsTransform() const { return Rivet::cmsTransform(beams()); }

    /// Get the beta factor vector for the Lorentz boost to the beam centre-of-mass
    Vector3 cmsBetaVec() const { return Rivet::cmsBetaVec(beams()); }

    /// Get the gamma factor vector for the Lorentz boost to the beam centre-of-mass
    Vector3 cmsGammaVec() const { return Rivet::cmsGammaVec(beams()); }

    /// @}


    /// @name Per-nucleon beam kinematics
    /// @{

    /// Get per-nucleon centre of mass energy, \f$ \sqrt{s}/(A_1 + A_2) \f$
    double asqrtS() const { return Rivet::asqrtS(beams()); }

    /// Get the Lorentz boost to the per-nucleon beam centre-of-mass
    Vector3 acmsBetaVec() const { return Rivet::acmsBetaVec(beams()); }

    /// Get the Lorentz boost to the per-nucleon beam centre-of-mass
    Vector3 acmsGammaVec() const { return Rivet::acmsGammaVec(beams()); }

    /// Get the Lorentz transform to the per-nucleon beam centre-of-mass
    LorentzTransform acmsTransform() const { return Rivet::acmsTransform(beams()); }

    /// @}


    /// Get the beam interaction primary vertex (PV) position
    FourVector pv() const;


    /// Project on to the Event
    virtual void project(const Event& e);


  protected:

    /// The beam particles in the current collision
    ParticlePair _theBeams;

  private:

    /// Compare with other projections -- it's always the same, since there are no params
    virtual CmpState compare(const Projection&) const { return CmpState::EQ; }

  };


}

#endif
