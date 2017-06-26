// -*- C++ -*-
#ifndef RIVET_SmearedParticles_HH
#define RIVET_SmearedParticles_HH

#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Tools/SmearingFunctions.hh"

namespace Rivet {


  // Recursive variadic template decoding
  namespace {
    // For handling an optional Cut as final element... too fiddly?
    // void toEffSmearFns(vector<ParticleEffSmearFn>&, Cut& c, const Cut& c2) {
    //   c = c2;
    // }
    template<typename T>
    void toEffSmearFns(vector<ParticleEffSmearFn>& v, const T& t) {
      v.push_back(ParticleEffSmearFn(t));
    }
    template<typename T, typename... Args>
    void toEffSmearFns(vector<ParticleEffSmearFn>& v, const T& first, Args... args) {
      v.push_back(ParticleEffSmearFn(first));
      toEffSmearFns(v, args...);
    }
  }


  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  class SmearedParticles : public ParticleFinder {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with const efficiency
    SmearedParticles(const ParticleFinder& pf,
                     double eff,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, ParticleEffFilter(eff), PARTICLE_SMEAR_IDENTITY, c)
    {    }

    /// @brief Constructor with an efficiency function arg
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleEffFn& effFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, effFn, PARTICLE_SMEAR_IDENTITY, c)
    {    }

    /// @brief Constructor with const efficiency and a smearing function
    SmearedParticles(const ParticleFinder& pf,
                     double eff, const ParticleSmearFn& smearFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, ParticleEffFilter(eff), smearFn, c)
    {    }

    /// @brief Constructor with efficiency and smearing function args
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleEffFn& effFn, const ParticleSmearFn& smearFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, c, effFn, smearFn)
    {    }

    // /// @brief Constructor with efficiency and smearing function args
    // SmearedParticles(const ParticleFinder& pf,
    //                  const ParticleEffFn& effFn, const ParticleSmearFn& smearFn,
    //                  const Cut& c=Cuts::open())
    //   : ParticleFinder(c),
    //     _effFn(effFn), _smearFn(smearFn)
    // {
    //   setName("SmearedParticles");
    //   addProjection(pf, "TruthParticles");
    // }

    /// @brief Constructor with efficiency and smearing function args
    SmearedParticles(const ParticleFinder& pf,
                     const vector<ParticleEffSmearFn>& effSmearFns,
                     const Cut& c=Cuts::open())
      : ParticleFinder(c),
        _detFns(effSmearFns)
    {
      setName("SmearedParticles");
      addProjection(pf, "TruthParticles");
    }

    /// @brief Constructor with efficiency and smearing function args
    SmearedParticles(const ParticleFinder& pf,
                     const initializer_list<ParticleEffSmearFn>& effSmearFns,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, vector<ParticleEffSmearFn>{effSmearFns}, c)
    {    }

    /// @brief Constructor with an arbitrary list of efficiency and smearing function args
    /// @todo Wouldn't it be nice if the Cut could go *after* the parameter pack? Oh well...
    template<typename... Args>
    SmearedParticles(const ParticleFinder& pf, const Cut& c, Args... effSmearFns)
      : ParticleFinder(c)
    {
      setName("SmearedParticles");
      addProjection(pf, "TruthParticles");

      toEffSmearFns(_detFns, effSmearFns...);
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SmearedParticles);

    //@}


    /// Compare to another SmearedParticles
    int compare(const Projection& p) const {
      /// @todo Need to implement the comparison operator on the *EffSmearFn structs
      return UNDEFINED;

      /// STD::FUNCTION VERSION
      // const SmearedParticles& other = dynamic_cast<const SmearedParticles&>(p);
      // if (get_address(_detFns[0]) == 0) return UNDEFINED;
      // MSG_TRACE("hashes = ");
      // for (size_t i = 0; i < _detFns.size(); ++i)
      //   MSG_TRACE( get_address(_detFns[i]) << "," << get_address(other._detFns[i]) << "; ");
      // Cmp<unsigned long> ret = mkPCmp(other, "TruthParticles");
      // for (size_t i = 0; i < _detFns.size(); ++i)
      //   ret = ret || cmp(get_address(_detFns[i]), get_address(other._detFns[i]));
      // return ret;

      /// OLD VERSION
      // if (get_address(_effFn) == 0) return UNDEFINED;
      // if (get_address(_smearFn) == 0) return UNDEFINED;
      // MSG_TRACE("Eff hashes = " << get_address(_effFn) << "," << get_address(other._effFn) << "; " <<
      //           "smear hashes = " << get_address(_smearFn) << "," << get_address(other._smearFn));
      // return mkPCmp(other, "TruthParticles") ||
      //   cmp(get_address(_effFn), get_address(other._effFn)) ||
      //   cmp(get_address(_smearFn), get_address(other._smearFn));
    }


    /// Perform the particle finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Particles& truthparticles = apply<ParticleFinder>(e, "TruthParticles").particlesByPt();
      _theParticles.clear(); _theParticles.reserve(truthparticles.size());
      for (const Particle& p : truthparticles) {
        Particle pdet = p;
        double peff = -1;
        bool keep = true;
        for (const ParticleEffSmearFn& fn : _detFns) {
          tie(pdet, peff) = fn(pdet); // smear & eff
          MSG_DEBUG("New det particle: pid=" << pdet.pid()
                    << ", mom=" << pdet.mom()/GeV << " GeV, "
                    << "pT=" << pdet.pT()/GeV << ", eta=" << pdet.eta()
                    << " : eff=" << 100*peff << "%");
          if (peff <= 0) { keep = false; break; } //< no need to roll expensive dice (and we deal with -ve probabilities, just in case)
          if (peff < 1 && rand01() > peff)  { keep = false; break; } //< roll dice (and deal with >1 probabilities, just in case)
        }
        if (keep) _theParticles.push_back(pdet);
      }
    }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _theParticles.clear(); }


  private:

    /// Stored efficiency & smearing functions
    vector<ParticleEffSmearFn> _detFns;

  };


}

#endif
