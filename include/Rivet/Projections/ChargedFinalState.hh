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
  /// @todo Can we make this work nicely with inheritance rather than composition?
  class ChargedFinalState : public FinalState {

  public:
    
    /// @name Constructors
    //@{
    ChargedFinalState(const FinalState& fsp) { 
      setName("ChargedFinalState");
      addProjection(fsp, "FS");
    }
    
    ChargedFinalState(double mineta = -MaxRapidity,
                      double maxeta =  MaxRapidity,
                      double minpt  =  0.0*GeV)
    { 
      setName("ChargedFinalState");
      addProjection(*new FinalState(mineta, maxeta, minpt), "FS");
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
