// -*- C++ -*-
#ifndef RIVET_FinalState_HH
#define RIVET_FinalState_HH

#include "Rivet/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet {


  /// Project out all final-state particles in an event.
  class FinalState : public Projection {    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    FinalState(double mineta = -MAXRAPIDITY,
               double maxeta =  MAXRAPIDITY,
               double minpt  =  0.0*GeV);
    
    /// A constructor which allows to specify multiple eta ranges
    /// and the min \f$ p_T \f$ (in GeV).
    FinalState(const vector<pair<double, double> >& etaRanges,
               double minpt = 0.0*GeV);

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new FinalState(*this);
    }

    //@}


    /// Get the final-state particles.
    virtual const ParticleVector& particles() const { return _theParticles; }

    /// Get the final-state particles, ordered by supplied sorting function object.
    template <typename F>
    const ParticleVector& particles(F sorter) const {
      std::sort(_theParticles.begin(), _theParticles.end(), sorter);
      return _theParticles;
    }

    /// Get the final-state particles, ordered by \f$ p_T \f$.
    const ParticleVector& particlesByPt() const {
      return particles(cmpParticleByPt);
    }

    /// Get the final-state particles, ordered by \f$ E \f$.
    const ParticleVector& particlesByE() const {
      return particles(cmpParticleByE);
    }

    /// Get the final-state particles, ordered by \f$ E_T \f$.
    const ParticleVector& particlesByEt() const {
      return particles(cmpParticleByEt);
    }


    /// Access the projected final-state particles.
    virtual const size_t size() const { return _theParticles.size(); }

    /// Is this final state empty?
    virtual const bool empty() const { return _theParticles.empty(); }


  public:

    typedef Particle entity_type;
    typedef ParticleVector collection_type; 

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
 
    /// The ranges allowed for pseudorapidity.
    vector<pair<double,double> > _etaRanges;
    
    /// The minimum allowed transverse momentum.
    double _ptmin;
    
    /// The final-state particles.
    mutable ParticleVector _theParticles;
    
  };

  
}

#endif
