// -*- C++ -*-
#ifndef RIVET_SmearedParticles_HH
#define RIVET_SmearedParticles_HH

#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Tools/SmearingFunctions.hh"

namespace Rivet {


  // Recursive variadic template arg decoding
  namespace {
    template<typename T>
    vector<ParticleEffSmearFn>& toEffSmearFns(vector<ParticleEffSmearFn>& v, const T& t) {
      v.push_back(ParticleEffSmearFn(t));
      return v;
    }
    template<typename T, typename... ARGS>
    vector<ParticleEffSmearFn>& toEffSmearFns(vector<ParticleEffSmearFn>& v, const T& first, ARGS... args) {
      v.push_back(ParticleEffSmearFn(first));
      toEffSmearFns(v, args...);
      return v;
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
      : SmearedParticles(pf, {eff}, c)
    {    }

    /// @brief Constructor with an efficiency function
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleEffFn& effFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {effFn}, c)
    {    }

    /// @brief Constructor with const efficiency followed by a smearing function
    SmearedParticles(const ParticleFinder& pf,
                     double eff, const ParticleSmearFn& smearFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {eff, smearFn}, c)
    {    }

    /// @brief Constructor with a smearing function followed by const efficiency
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleSmearFn& smearFn, double eff,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {smearFn, eff}, c)
    {    }

    /// @brief Constructor with an efficiency function followed by a smearing function
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleEffFn& effFn, const ParticleSmearFn& smearFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {effFn, smearFn}, c)
    {    }

    /// @brief Constructor with a smearing function followed by an efficiency function
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleSmearFn& smearFn, const ParticleEffFn& effFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {smearFn, effFn}, c)
    {    }

    /// @brief Constructor with an ordered list of efficiency and/or smearing functions
    SmearedParticles(const ParticleFinder& pf,
                     const vector<ParticleEffSmearFn>& effSmearFns,
                     const Cut& c=Cuts::open())
      : ParticleFinder(c),
        _detFns(effSmearFns)
    {
      setName("SmearedParticles");
      addProjection(pf, "TruthParticles");
    }

    /// @brief Constructor with an ordered list of efficiency and/or smearing functions
    SmearedParticles(const ParticleFinder& pf,
                     const initializer_list<ParticleEffSmearFn>& effSmearFns,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, vector<ParticleEffSmearFn>{effSmearFns}, c)
    {    }

    /// @brief Constructor with a variadic ordered list of efficiency and smearing function args
    /// @note The Cut must be provided *before* the eff/smearing functions
    /// @todo Wouldn't it be nice if the Cut could also go *after* the parameter pack?
    template<typename... ARGS>
    SmearedParticles(const ParticleFinder& pf, const Cut& c, ARGS... effSmearFns)
      : SmearedParticles(pf, toEffSmearFns(_detFns, effSmearFns...), c)
    {    }


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
