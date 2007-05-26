// -*- C++ -*-
#ifndef RIVET_TotalVisibleMomentum_H
#define RIVET_TotalVisibleMomentum_H

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/RivetCLHEP.fhh"
#include "Rivet/Event.hh"


namespace Rivet {

  /// Project out the total visible energy vector, allowing missing ET etc to be calculated.
  class TotalVisibleMomentum: public Projection {
    
  public:
    
    /// Constructor. The provided FinalState projection must live throughout the run.
    inline TotalVisibleMomentum(VetoedFinalState& vfsp)
      : vfsproj(&vfsp)
    { 
      addProjection(vfsp);
    }
    
  public:
    /// Return the name of the projection
    inline string getName() const {
      return "TotalVisibleMomentum";
    }

    /// The projected four-momentum vector
    inline CLHEP::LorentzVector& getMomentum() { return _momentum; }

    /// The projected four-momentum vector
    inline const CLHEP::LorentzVector& getMomentum() const { return _momentum; }

    
  protected:
    
    /// Apply the projection to the event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection & p) const;
        
  private:
        
    /// The FinalState projection used by this projection
    VetoedFinalState * vfsproj;

    /// The total visible momentum
    CLHEP::LorentzVector _momentum;
    
  private:
    
    /// Hide the assignment operator.
    TotalVisibleMomentum & operator=(const TotalVisibleMomentum &);
    
  };
  
}


#endif
