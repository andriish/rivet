// -*- C++ -*-
#ifndef RIVET_FinalState_HH
#define RIVET_FinalState_HH

#include "Rivet/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Cuts.hh"


namespace Rivet {


  /// @brief Project out all final-state particles in an event.
  /// Probably the most important projection in Rivet!
  class FinalState : public Projection {
  public:

    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    FinalState(double mineta,
               double maxeta,
               double minpt = 0.0);

    /// Testing construction using Cuts object
    FinalState(Cut c = Cuts::open());

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new FinalState(*this);
    }

    //@}


    /// Get the final-state particles.
    virtual const Particles& particles() const { return _theParticles; }

    /// Get the final-state particles, ordered by supplied sorting function object.
    template <typename F>
    const Particles& particles(F sorter) const {
      std::sort(_theParticles.begin(), _theParticles.end(), sorter);
      return _theParticles;
    }

    /// Get the final-state particles, ordered by decreasing \f$ p_T \f$.
    const Particles& particlesByPt() const {
      return particles(cmpParticleByPt);
    }

    /// Get the final-state particles, ordered by decreasing \f$ p \f$.
    const Particles& particlesByP() const {
      return particles(cmpParticleByP);
    }

    /// Get the final-state particles, ordered by decreasing \f$ E \f$.
    const Particles& particlesByE() const {
      return particles(cmpParticleByE);
    }

    /// Get the final-state particles, ordered by decreasing \f$ E_T \f$.
    const Particles& particlesByEt() const {
      return particles(cmpParticleByEt);
    }

    /// Get the final-state particles, ordered by increasing \f$ \eta \f$.
    const Particles& particlesByEta() const {
      return particles(cmpParticleByAscPseudorapidity);
    }

    /// Get the final-state particles, ordered by increasing \f$ |\eta| \f$.
    const Particles& particlesByModEta() const {
      return particles(cmpParticleByAscAbsPseudorapidity);
    }

    /// Get the final-state particles, ordered by increasing \f$ y \f$.
    const Particles& particlesByRapidity() const {
      return particles(cmpParticleByAscRapidity);
    }

    /// Get the final-state particles, ordered by increasing \f$ |y| \f$.
    const Particles& particlesByModRapidity() const {
      return particles(cmpParticleByAscAbsRapidity);
    }

    /// Access the projected final-state particles.
    virtual size_t size() const { return _theParticles.size(); }

    /// Is this final state empty?
    virtual bool empty() const { return _theParticles.empty(); }
    /// @deprecated Is this final state empty?
    virtual bool isEmpty() const { return _theParticles.empty(); }

    /// Minimum-\f$ p_\perp \f$ requirement.
    //virtual double ptMin() const { return _ptmin; }


  public:

    typedef Particle entity_type;
    typedef Particles collection_type;

    /// Template-usable interface common to JetAlg.
    const collection_type& entities() const {
      return particles();
    }


  protected:

    /// Apply the projection to the event.
    virtual void project(const Event& e);

    /// Compare projections.
    virtual int compare(const Projection& p) const;

    /// Decide if a particle is to be accepted or not.
    bool accept(const Particle& p) const;


  protected:
    /// The applicable cuts
    Cut _cuts;

    /// The final-state particles.
    mutable Particles _theParticles;

  };


}

#endif
