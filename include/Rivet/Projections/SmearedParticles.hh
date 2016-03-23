// -*- C++ -*-
#ifndef RIVET_SmearedParticles_HH
#define RIVET_SmearedParticles_HH

#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include "boost/function.hpp"

namespace Rivet {


  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  class SmearedParticles : public ParticleFinder {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with efficiency and smearing function args
    template <typename P2DFN, typename P2PFN>
    SmearedParticles(const ParticleFinder& pf,
                     const P2DFN& effFn, const P2PFN& smearFn=PARTICLE_SMEAR_IDENTITY,
                     const Cut& c=Cuts::open())
      : ParticleFinder(c),
        _effFnPtr(reinterpret_cast<size_t>(effFn)), _smearFnPtr(reinterpret_cast<size_t>(smearFn)),
        _effFn(effFn), _smearFn(smearFn)
    {
      setName("SmearedParticles");
      addProjection(pf, "TruthParticles");
    }


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new SmearedParticles(*this);
    }

    //@}


    /// Compare to another SmearedParticles
    int compare(const Projection& p) const {
      const SmearedParticles& other = dynamic_cast<const SmearedParticles&>(p);
      return cmp(_effFnPtr, other._effFnPtr) || cmp(_smearFnPtr, other._smearFnPtr);
    }


    /// Perform the particle finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Particles& truthparticles = applyProjection<ParticleFinder>(e, "TruthParticles").particlesByPt();
      _theParticles.clear(); _theParticles.reserve(truthparticles.size());
      foreach (const Particle& p, truthparticles) {
        const double peff = (_effFn) ? _effFn(p) : 1;
        MSG_DEBUG("Efficiency of particle " << p.mom() << " = " << 100*peff << "%");
        if (peff == 0) continue; //< no need to roll expensive dice
        if (peff == 1 || peff < rand01()) {
          _theParticles.push_back(_smearFn ? _smearFn(p) : p); //< smearing
        }
      }
    }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _theParticles.clear(); }


  private:

    // Particles _recoparticles;
    size_t _effFnPtr, _smearFnPtr;
    boost::function<double(const Particle&)> _effFn;
    boost::function<Particle(const Particle&)> _smearFn;

  };


}

#endif
