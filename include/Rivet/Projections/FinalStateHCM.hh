// -*- C++ -*-
#ifndef RIVET_FinalStateHCM_H
#define RIVET_FinalStateHCM_H

#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"


namespace Rivet {

  /// Project all final state particles (except the scattered lepton)
  /// boosted to the hadronic center of mass system.
  class FinalStateHCM: public FinalState {

  public:
    
    /// The default constructor. Must specify DISLepton, DISKinematics
    /// and FinalState projection objects which are assumed to live
    /// throughout the run.
    inline FinalStateHCM(DISLepton& leptonp, DISKinematics& kinematicsp, FinalState& fsp)
      : _lepton(&leptonp), _kinematics(&kinematicsp), _fsproj(&fsp) 
    { 
      addProjection(leptonp);
      addProjection(kinematicsp);
      addProjection(fsp);
    }
    

  public:
    /// Return the name of the projection
    inline string getName() const {
      return "FinalStateHCM";
    }
    
  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
    
  public:
    
    /// Access the projected final-state particles.
    inline const ParticleVector& particles() const { return _theParticles; }
    
  private:
    
    /// The projector for the DIS lepton.
    DISLepton* _lepton;
    
    /// The projector for the DIS kinematics.
    DISKinematics* _kinematics;

    /// The projector for the full final state.
    FinalState* _fsproj;
    
    /// The final-state particles.
    ParticleVector _theParticles;
    
  private:
    
    /// Hide the assignment operator.
    FinalStateHCM& operator=(const FinalStateHCM&);
    
  };
  
}


#endif
