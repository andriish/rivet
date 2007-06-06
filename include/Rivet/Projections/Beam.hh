// -*- C++ -*-
#ifndef RIVET_Beam_H
#define RIVET_Beam_H

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"

namespace Rivet {

  /// Project out the incoming beams
  class Beam : public Projection {

  public:

    /// The default constructor.
    inline Beam() { }

  public:
      /// Return the name of the projection
      inline string getName() const {
        return "Beam";
      }

  public:

    /// The pair of beam particles in the current collision in GenEvent 
    inline const ParticlePair& operator()() const {
      return theBeams;
    }

  protected:

    /// Project on to the Event
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual int compare(const Projection & p) const;


  private:

    /// The beam particles in the current collision in GenEvent 
    ParticlePair theBeams;

    /// Hiding the assignment operator.
    Beam & operator=(const Beam &);

  };

}


#endif /* RIVET_Beam_H */
