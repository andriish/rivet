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
    ChargedFinalState(FinalState& fsp)
      : FinalState(fsp)
    { }
    
    ChargedFinalState(double mineta = -MaxRapidity,
                      double maxeta = MaxRapidity,
                      double minpt = 0.0)
      : FinalState(mineta, maxeta, minpt)
    { }

    ~ChargedFinalState() {
      getLog() << Log::TRACE << "Destroying " << getName() << " at " << this << endl;
    }


    /// Return the name of the projection.
    string getName() const {
      return "ChargedFinalState";
    }

  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;

  private:

    /// @todo For now, we'll hold a constituent FinalState...
    FinalState _fsproj;
    
  };

  
}


#endif
