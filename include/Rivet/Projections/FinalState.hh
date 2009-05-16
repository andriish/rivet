// -*- C++ -*-
#ifndef RIVET_FinalState_HH
#define RIVET_FinalState_HH

#include "Rivet/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet {

  /// Project out all final-state particles in an event.
  class FinalState: public Projection {
    
  public:
    
    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    FinalState(double mineta = -MaxRapidity,
               double maxeta =  MaxRapidity,
               double minpt  =  0.0*GeV);
    
    /// A constructor which allows to specify multiple eta ranges
    /// and the min \f$ p_T \f$ (in GeV).
    FinalState(const std::vector<std::pair<double, double> >& etaRanges,
               double minpt  =  0.0*GeV);

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new FinalState(*this);
    }

    //@}


    /// Access the projected final-state particles.
    virtual const ParticleVector& particles() const { return _theParticles; }

    /// Access the projected final-state particles.
    virtual const size_t size() const { return _theParticles.size(); }

    /// Is this final state empty?
    virtual const bool isEmpty() const { return _theParticles.empty(); }


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
    std::vector<std::pair<double,double> > _etaRanges;
    
    /// The minimum allowed transverse momentum.
    double _ptmin;
    
    /// The final-state particles.
    ParticleVector _theParticles;
    
  };
  
}


#endif
