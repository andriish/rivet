// -*- C++ -*-
#ifndef RIVET_Beam_HH
#define RIVET_Beam_HH

#include "Rivet/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"


namespace Rivet {

  /// Project out the incoming beams
  class Beam : public Projection {
    
  public:
    
    /// The default constructor.
    Beam() { }

    /// The default destructor.
    ~Beam() { }

    /// Return the name of the projection
    string getName() const {
      return "Beam";
    }
    
    /// The pair of beam particles in the current collision in GenEvent 
    const ParticlePair& getBeams() const {
      return _theBeams;
    }


  protected:
    /// Project on to the Event
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual int compare(const Projection& p) const {
      return 0;
    }


  private:
    /// The beam particles in the current collision in GenEvent 
    ParticlePair _theBeams;

  };

}

#endif
