// -*- C++ -*-
#ifndef RIVET_HadronicFinalState_HH
#define RIVET_HadronicFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {

  /// Project only charged final state particles.
  class HadronicFinalState : public FinalState {

  public:
    
    /// Constructor: the supplied FinalState projection is assumed to live through the run.
    HadronicFinalState(FinalState& fsp) 
      : FinalState(fsp)
    { 
      setName("HadronicFinalState");
      addProjection(fsp, "FS");
    }
    
    HadronicFinalState(double mineta = -MaxRapidity,
                       double maxeta = MaxRapidity,
                       double minpt = 0.0*GeV)
      : FinalState(mineta, maxeta, minpt)
    {
      setName("HadronicFinalState");
      addProjection(FinalState(mineta, maxeta, minpt), "FS");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new HadronicFinalState(*this);
    }

  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
    
  };

  
}


#endif
