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
    Beam() { 
      setName("Beam");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new Beam(*this);
    }


  public:

    /// The pair of beam particles in the current collision.
    const ParticlePair& beams() const {
      return _theBeams;
    }

    /// The pair of beam particle PDG codes in the current collision.
    const pair<long,long> beamIDs() const {
      return make_pair(beams().first.pdgId(), 
                       beams().second.pdgId());
    }

    /// Get centre of mass energy, \f$ \sqrt{s} \f$.
    const double sqrtS() const;


  protected:
    /// Project on to the Event
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual int compare(const Projection& p) const {
      return PCmp::EQUIVALENT;
    }


  private:
    /// The beam particles in the current collision in GenEvent 
    ParticlePair _theBeams;

  };

}

#endif
