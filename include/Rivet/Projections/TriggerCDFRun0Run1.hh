// -*- C++ -*-
#ifndef RIVET_TriggerCDFRun0Run1_HH
#define RIVET_TriggerCDFRun0Run1_HH

#include "Rivet/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// Project out the incoming beams
  class TriggerCDFRun0Run1 : public Projection {
  public:
    
    /// Default constructor.
    TriggerCDFRun0Run1() { 
      setName("TriggerCDFRun0Run1");

      addProjection(ChargedFinalState(-5.9, 5.9), "CFS");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new TriggerCDFRun0Run1(*this);
    }


  public:

    /// The trigger result
    const bool minBiasDecision() const {
      return _decision_mb;
    }

    /// Project on to the Event
    void project(const Event& evt);


  protected:

    /// Compare with other projections.
    virtual int compare(const Projection& p) const {
      return EQUIVALENT;
    }


  private:

    /// The min bias trigger decision
    bool _decision_mb;

  };


}

#endif
