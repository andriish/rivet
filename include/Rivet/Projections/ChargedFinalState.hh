// -*- C++ -*-
#ifndef RIVET_ChargedFinalState_HH
#define RIVET_ChargedFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// Project only charged final state particles.
  class ChargedFinalState : public FinalState {

  public:
    
    /// @name Constructors
    //@{
    ChargedFinalState(const FinalState& fsp) { 
      setName("ChargedFinalState");
      addProjection(fsp, "FS");
    }
    
    ChargedFinalState(double mineta = -MAXRAPIDITY,
                      double maxeta =  MAXRAPIDITY,
                      double minpt  =  0.0*GeV)
    { 
      setName("ChargedFinalState");
      addProjection(FinalState(mineta, maxeta, minpt), "FS");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new ChargedFinalState(*this);
    }
    //@}


  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
  };

  
}


#endif
