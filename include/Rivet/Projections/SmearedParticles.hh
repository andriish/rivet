// -*- C++ -*-
#ifndef RIVET_SmearedParticles_HH
#define RIVET_SmearedParticles_HH

#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include <functional>

namespace Rivet {


  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  class SmearedParticles : public ParticleFinder {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with efficiency and smearing function args
    template <typename P2DFN>
    SmearedParticles(const ParticleFinder& pf,
                     const P2DFN& effFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, effFn, PARTICLE_SMEAR_IDENTITY, c)
    {    }


    /// @brief Constructor with efficiency and smearing function args
    template <typename P2DFN, typename P2PFN>
    SmearedParticles(const ParticleFinder& pf,
                     const P2DFN& effFn, const P2PFN& smearFn,
                     const Cut& c=Cuts::open())
      : ParticleFinder(c),
        _effFn(effFn), _smearFn(smearFn)
    {
      setName("SmearedParticles");
      addProjection(pf, "TruthParticles");
    }


    /// Clone on the heap.
    virtual unique_ptr<Projection> clone() const {
      return unique_ptr<Projection>(new SmearedParticles(*this));
    }

    //@}


    /// Compare to another SmearedParticles
    int compare(const Projection& p) const {
      const SmearedParticles& other = dynamic_cast<const SmearedParticles&>(p);
      if (_mkhash(_effFn) == 0) return UNDEFINED;
      if (_mkhash(_smearFn) == 0) return UNDEFINED;
      return cmp(_mkhash(_effFn), _mkhash(other._effFn)) || cmp(_mkhash(_smearFn), _mkhash(other._smearFn));
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

    /// Make a hash integer from the provided wrapped Particle -> double function
    size_t _mkhash(const std::function<double(const Particle&)>& fn) const {
      const size_t rtn = reinterpret_cast<size_t>(fn.target<double(*)(const Particle&)>());
      MSG_TRACE("P2D hash = " << rtn);
      return rtn;
    }

    /// Make a hash integer from the provided wrapped Particle -> Particle function
    size_t _mkhash(const std::function<Particle(const Particle&)>& fn) const {
      const size_t rtn = reinterpret_cast<size_t>(fn.target<Particle(*)(const Particle&)>());
      MSG_TRACE("P2P hash = " << rtn);
      return rtn;
    }

    /// Stored efficiency function
    std::function<double(const Particle&)> _effFn;
    /// Stored smearing function
    std::function<Particle(const Particle&)> _smearFn;

  };


}

#endif
