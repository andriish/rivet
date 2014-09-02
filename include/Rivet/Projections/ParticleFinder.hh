// -*- C++ -*-
#ifndef RIVET_ParticleFinder_HH
#define RIVET_ParticleFinder_HH

#include "Rivet/Projection.hh"
#include "Rivet/Cuts.hh"

namespace Rivet {


  /// @brief Base class for projections which return subsets of an event's particles
  class ParticleFinder : public Projection {
  public:

    /// @name Object lifetime management
    //@{

    /// Construction using Cuts object
    ParticleFinder(Cut c = Cuts::open()) : _cuts(c), _theParticles() {}

    /// Virtual destructor for inheritance
    virtual ~ParticleFinder() {}

    /// Clone on the heap.
    virtual const Projection* clone() const = 0;

    //@}


    /// @name Particle accessors
    //@{

    /// Access the projected final-state particles.
    virtual size_t size() const { return _theParticles.size(); }

    /// Is this final state empty?
    virtual bool empty() const { return _theParticles.empty(); }
    /// @deprecated Is this final state empty?
    virtual bool isEmpty() const { return _theParticles.empty(); }

    /// Get the final-state particles in no particular order, with no cuts.
    virtual const Particles& particles() const { return _theParticles; }

    /// @brief Get the final-state particles, with optional cuts.
    /// @note Returns a copy rather than a reference, due to cuts
    //virtual Particles particles(Cut c) const {  }

    /// Get the final-state particles, ordered by supplied sorting function object.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @todo Update to use Cuts
    template <typename F>
    Particles particles(Cut c, F sorter) const {
      Particles result;
      result.reserve(size());
      if ( c == Cuts::open() )
      	result.assign(_theParticles.begin(), _theParticles.end());
      else 
	foreach (const Particle& p, _theParticles)
      	  if ( c->accept(p) ) result.push_back(p);
      std::sort(result.begin(), result.end(), sorter);
      return result;
    }

    /// Get the final-state particles, ordered by decreasing \f$ p_T \f$.
    /// @todo Update to use Cuts
    Particles particlesByPt(Cut c = Cuts::open()) const {
      return particles(c, cmpMomByPt);
    }

    /// Get the final-state particles, ordered by decreasing \f$ p \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    Particles particlesByP(Cut c = Cuts::open()) const {
      return particles(c, cmpMomByP);
    }

    /// Get the final-state particles, ordered by decreasing \f$ E \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    Particles particlesByE(Cut c = Cuts::open()) const {
      return particles(c, cmpMomByE);
    }

    /// Get the final-state particles, ordered by decreasing \f$ E_T \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    Particles particlesByEt(Cut c = Cuts::open()) const {
      return particles(c, cmpMomByEt);
    }

    /// Get the final-state particles, ordered by increasing \f$ \eta \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    Particles particlesByEta(Cut c = Cuts::open()) const {
      return particles(c, cmpMomByAscPseudorapidity);
    }

    /// Get the final-state particles, ordered by increasing \f$ |\eta| \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    Particles particlesByModEta(Cut c = Cuts::open()) const {
      return particles(c, cmpMomByAscAbsPseudorapidity);
    }

    /// Get the final-state particles, ordered by increasing \f$ y \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    Particles particlesByRapidity(Cut c = Cuts::open()) const {
      return particles(c, cmpMomByAscRapidity);
    }

    /// Get the final-state particles, ordered by increasing \f$ |y| \f$.
    /// @todo Update to use Cuts
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    Particles particlesByModRapidity(Cut c = Cuts::open()) const {
      return particles(c, cmpMomByAscAbsRapidity);
    }

    //@}


    /// Minimum-\f$ p_\perp \f$ requirement.
    /// @todo Replace with cuts() accessor
    ///virtual Cut cuts() const { return _cuts; }


    /// @name JetAlg compatibility
    //@{

    typedef Particle entity_type;
    typedef Particles collection_type;

    /// Template-usable interface common to JetAlg.
    const collection_type& entities() const {
      return particles();
    }

    //@}


  protected:

    /// Apply the projection to the event.
    virtual void project(const Event& e) = 0;

    /// Compare projections.
    virtual int compare(const Projection& p) const;

  protected:
    /// The applicable cuts
    Cut _cuts;

    /// The final-state particles.
    Particles _theParticles;

  };


}

#endif
