// -*- C++ -*-
#ifndef RIVET_SmearedParticles_HH
#define RIVET_SmearedParticles_HH

#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include <functional>

namespace {
  /// Make a hash integer from an std::function
  template<typename T, typename... U>
  inline size_t getaddr(std::function<T(U...)> f) {
    typedef T(fnType)(U...);
    fnType ** fnPointer = f.template target<fnType*>();
    return (fnPointer != nullptr) ? reinterpret_cast<size_t>(*fnPointer) : 0;
  }
}

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
      if (getaddr(_effFn) == 0) return UNDEFINED;
      if (getaddr(_smearFn) == 0) return UNDEFINED;
      MSG_TRACE("Eff hashes = " << getaddr(_effFn) << "," << getaddr(other._effFn) << "; " <<
                "smear hashes = " << getaddr(_smearFn) << "," << getaddr(other._smearFn));
      return cmp(getaddr(_effFn), getaddr(other._effFn)) || cmp(getaddr(_smearFn), getaddr(other._smearFn));
    }


    /// Perform the particle finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Particles& truthparticles = applyProjection<ParticleFinder>(e, "TruthParticles").particlesByPt();
      _theParticles.clear(); _theParticles.reserve(truthparticles.size());
      for (const Particle& p : truthparticles) {
        const double peff = (_effFn) ? _effFn(p) : 1;
        MSG_DEBUG("Efficiency of particle with pid=" << p.pid()
                  << ", mom=" << p.mom()/GeV << " GeV, "
                  << "pT=" << p.pT()/GeV << ", eta=" << p.eta()
                  << " : " << 100*peff << "%");
        if (peff == 0) continue; //< no need to roll expensive dice
        if (peff == 1 || peff < rand01()) {
          _theParticles.push_back(_smearFn ? _smearFn(p) : p); //< smearing
        }
      }
    }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _theParticles.clear(); }


  private:

    /// Stored efficiency function
    std::function<double(const Particle&)> _effFn;

    /// Stored smearing function
    std::function<Particle(const Particle&)> _smearFn;

  };


}

#endif
