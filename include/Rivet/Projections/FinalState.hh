// -*- C++ -*-
#ifndef RIVET_FinalState_H
#define RIVET_FinalState_H

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet {

  /// Project out all final-state particles in an event.
  class FinalState: public Projection {
    
  public:
    
    /// @name Standard constructors and destructors. */
    //@{
    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    inline FinalState(double mineta = -MaxRapidity,
                      double maxeta = MaxRapidity,
                      double minpt = 0.0)
      : _etamin(mineta), _etamax(maxeta), _ptmin(minpt) 
    { 
      addCut("eta", MORE_EQ, mineta);
      addCut("eta", LESS_EQ, maxeta);
      addCut("pT",  MORE_EQ, minpt);
    }
    
  public:
    /// Return the name of the projection
    inline string getName() const {
      return "FinalState";
    }
    
  protected:
    
    /// Apply the projection to the event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
    
  public:
    
    /// Access the projected final-state particles.
    inline const ParticleVector& particles() const { return _theParticles; }
    
  private:
    
    /// The minimum allowed pseudo-rapidity.
    double _etamin;
    
    /// The maximum allowed pseudo-rapidity.
    double _etamax;
    
    /// The minimum allowed transverse momentum.
    double _ptmin;
    
    /// The final-state particles.
    ParticleVector _theParticles;
    
  private:
    
    /// Hide the assignment operator.
    FinalState& operator=(const FinalState&);
    
  };
  
}


#endif
