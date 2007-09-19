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
    
    /// Constructor: the supplied FinalState projection is assumed to live through the run.
    inline ChargedFinalState(FinalState& fsp)
      : _fsproj(fsp)
    {
      addProjection(_fsproj);
    }
    
    /// Return the name of the projection.
    inline string getName() const {
      return "ChargedFinalState";
    }

  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
    
  private:
    
    /// The projector for the full final state.
    FinalState _fsproj;
    
  };

  
}


#endif
