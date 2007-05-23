// -*- C++ -*-
#ifndef RIVET_VetoedFinalState_H
#define RIVET_VetoedFinalState_H

#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"


namespace Rivet {

  /// Project all final state particles except for those listed by PDG code.
  class VetoedFinalState : public FinalState {

  public:
    
    /// The default constructor. Must specify a FinalState projection 
    /// object which is assumed to live through the run.
    inline VetoedFinalState(FinalState& fsp, vector<long> vetocodes)
      : _fsproj(&fsp), _vetoCodes(vetocodes)
    {
      addProjection(fsp);
    }
    

  public:
    /// Return the name of the projection
    inline string getName() const {
      return "VetoedFinalState";
    }

  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
    
  public:
    
    /// Access the projected final-state particles.
    inline const ParticleVector& particles() const { return _theParticles; }

    /// Get the list of particle IDs to veto.
    inline const vector<long>& getVetoIds() const {
      return _vetoCodes;
    }

    /// Add a particle ID to veto.
    inline VetoedFinalState& addVetoId(const long id) {
      _vetoCodes.push_back(id);
      return *this;
    }

    /// Set the list of particle IDs to veto.
    inline VetoedFinalState& setVetoIds(const vector<long>& ids) {
      _vetoCodes = ids;
      return *this;
    }

    /// Clear the list of particle IDs to veto.
    inline VetoedFinalState& clearVetoIds() {
      _vetoCodes.clear();
      return *this;
    }
    
  private:
    
    /// The projector for the full final state.
    FinalState* _fsproj;
    
    /// The final-state particles.
    ParticleVector _theParticles;

    /// The final-state particles.
    vector<long> _vetoCodes;

    
  private:
    
    /// Hide the assignment operator.
    VetoedFinalState & operator=(const VetoedFinalState&);
    
  };
  
}


#endif
