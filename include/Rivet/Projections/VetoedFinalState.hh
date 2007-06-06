// -*- C++ -*-
#ifndef RIVET_VetoedFinalState_H
#define RIVET_VetoedFinalState_H
#
#include "Rivet/Tools/Logging.hh"
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
    inline VetoedFinalState(FinalState& fsp)
      : _fsproj(&fsp)
    {
      addProjection(fsp);
    }

    /// You can add a map of ID plus a 2D vector containing ptmin,
    /// ptmax defined the range of particles to be vetoed. 
    /// zero is the same as infinity for ptmax.
    inline VetoedFinalState(FinalState& fsp, const map<long,vector<double> >& vetocodes)
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

    /// Get the list of particle IDs and pt ranges to veto.
    inline const map<long,vector<double> >& getVetoIds() const {
      return _vetoCodes;
    }
  
    /// Add a particle ID and pt range to veto. ptmax=0.0 would give the same effect as infinity.
    inline VetoedFinalState& addVetoId(const long id, const double ptmin, const double ptmax) {
      vector<double> range; 
      range.push_back(ptmin);
      range.push_back(ptmax);
      _vetoCodes.insert(make_pair( id, range ));
      return *this;
    }

    /// Add a particle ID to veto (all pT range will be vetoed).
    inline VetoedFinalState& addVetoId(const long id) {
      vector<double> range; 
      _vetoCodes.insert(make_pair( id, range ));
      return *this;
    }

    /// Set the list of particle IDs and pt ranges to veto.
    inline VetoedFinalState& setVetoIds(const map<long,vector<double> >& ids) {
      _vetoCodes = ids;
      return *this;
    }

    /// Clear the list of particle IDs and ranges to veto.
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
    map<long,vector<double> > _vetoCodes;

    
  private:
    
    /// Hide the assignment operator.
    VetoedFinalState & operator=(const VetoedFinalState&);
    
  };
  
}


#endif
