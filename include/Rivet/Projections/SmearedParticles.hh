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

    /// @brief Constructor with vector of efficiency and smearing functions
    SmearedParticles(const ParticleFinder& pf,
                     const std::vector<std::function<pair<double,Particle>(const pair<double,Particle>&)> > detFn,
                     const Cut& c=Cuts::open())
      : ParticleFinder(c),
        _detFn(detFn)
    {
      setName("SmearedParticles");
      addProjection(pf, "TruthParticles");
    }

    template <typename DP2DPFN>
    /// @brief Constructor with efficiency or smearing function.
    /// Accepts lambdas, e.g. [](const std::pair<double,Particle>& p) { return make_pair(p.first,p.second) };
    SmearedParticles(const ParticleFinder& pf,
                     const DP2DPFN detFn,
                     const Cut& c=Cuts::open())
      : ParticleFinder(c)
    {
      _detFn.resize(1);
      _detFn[0] = detFn;
      setName("SmearedParticles");
      addProjection(pf, "TruthParticles");
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SmearedParticles);

    //@}


    /// Compare to another SmearedParticles
    int compare(const Projection& p) const {
      const SmearedParticles& other = dynamic_cast<const SmearedParticles&>(p);
      if (get_address(_detFn[0]) == 0) return UNDEFINED;
      MSG_TRACE("hashes = ");
      for(size_t i = 0; i < _detFn.size(); ++i) MSG_TRACE( get_address(_detFn[i]) << "," << get_address(other._detFn[i]) << "; ");
      Cmp<unsigned long> ret = mkPCmp(other, "TruthParticles");
      for(size_t i = 0; i < _detFn.size(); ++i) ret = ret || cmp(get_address(_detFn[i]), get_address(other._detFn[i]));
      return ret;
    }


    /// Perform the particle finding & smearing calculation
    void project(const Event& e) {
      // Filter the particles using given functions
      const Particles& truthparticles = apply<ParticleFinder>(e, "TruthParticles").particlesByPt();
      _theParticles.clear(); _theParticles.reserve(truthparticles.size());
      for (const Particle& p : truthparticles) {
        pair<double, Particle> tmp = make_pair(1., p);
        for(std::function<pair<double,Particle>(const pair<double,Particle>&)> fn : _detFn) tmp = fn(tmp);
        if(tmp.first <= 0) continue;
        if ( !(tmp.first < 1 && rand01() > tmp.first) ) _theParticles.push_back(tmp.second);
      }
    }

    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _theParticles.clear(); }


  private:

    /// Stored efficiency and smearing functions
    std::vector<std::function<pair<double,Particle>(const pair<double,Particle>&)> > _detFn;

  };


}

#endif
