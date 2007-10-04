// -*- C++ -*-
#ifndef RIVET_TotalVisibleMomentum_HH
#define RIVET_TotalVisibleMomentum_HH

#include "Rivet/Rivet.hh"
#include "Rivet/RivetCLHEP.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet {

  /// Project out the total visible energy vector, allowing missing 
  /// \f$ E_T \f$ etc. to be calculated.
  class TotalVisibleMomentum: public Projection {
    
  public:
    
    /// Constructor. The provided FinalState projection must live throughout the run.
    inline TotalVisibleMomentum(FinalState& fsp)
      : _fsproj(fsp)
    { 
      addProjection(fsp);
    }

    
  public:
    /// Return the name of the projection
    inline string getName() const {
      return "TotalVisibleMomentum";
    }

    /// The projected four-momentum vector
    inline LorentzVector& getMomentum() { return _momentum; }

    /// The projected four-momentum vector
    inline const LorentzVector& getMomentum() const { return _momentum; }

    /// The projected Scalar Transverse Momentum
    inline const double getSET() const { return _set; }

    
  protected:
    
    /// Apply the projection to the event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
        
  private:
    
    /// The FinalState projection used by this projection
    FinalState _fsproj;
    
    /// The total visible momentum
    LorentzVector _momentum;
    
    /// Scalar Transverse Energy
    double _set;
    
  };
  
}


#endif
