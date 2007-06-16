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
    
    /// Typedef for a vetoing entry.
    typedef map<long, pair<double, double> > VetoDetails;

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
    inline VetoedFinalState(FinalState& fsp, const VetoDetails& vetocodes)
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
    
    /// Get the list of particle IDs and \f$ p_T \f$ ranges to veto.
    inline const VetoDetails& getVetoDetails() const {
      return _vetoCodes;
    }
  
    /// Add a particle ID and \f$ p_T \f$ range to veto. Particles with \f$ p_T \f$ 
    /// IN the given range will be rejected.
    inline VetoedFinalState& addVetoDetail(const long id, const double ptmin, const double ptmax) {
      pair<double, double> ptrange; 
      ptrange.first = ptmin;
      ptrange.second = ptmax;
      _vetoCodes.insert(make_pair(id, ptrange));
      return *this;
    }

    /// Add a particle/antiparticle pair to veto in a given \f$ p_T \f$ range. Given a single ID, both
    /// the particle and its conjugate antiparticle will be rejected if their \f$ p_T \f$ is IN the given range.
    inline VetoedFinalState& addVetoPairDetail(const long id, const double ptmin, const double ptmax) {
      addVetoPairDetail(id,  ptmin, ptmax);
      addVetoPairDetail(-id, ptmin, ptmax);
      return *this;
    }

    /// Add a particle/antiparticle pair to veto. Given a single ID, both the particle and its corresponding 
    /// antiparticle (for all \f$ p_T \f$ values) will be vetoed.
    inline VetoedFinalState& addVetoPairId(const long id) {
      addVetoId(id);
      addVetoId(-id);
      return *this;
    }

    /// Add a particle ID to veto (all \f$ p_T \f$ range will be vetoed).
    inline VetoedFinalState& addVetoId(const long id) {
      pair<double, double> ptrange;
      ptrange.first = 0.0;
      ptrange.second = numeric_limits<double>::max();
      _vetoCodes.insert(make_pair(id, ptrange));
      return *this;
    }

    /// Set the list of particle IDs and \f$ p_T \f$ ranges to veto.
    inline VetoedFinalState& setVetoDetails(const VetoDetails& ids) {
      _vetoCodes = ids;
      return *this;
    }

    /// Clear the list of particle IDs and ranges to veto.
    inline VetoedFinalState& clearVetoDetails() {
      _vetoCodes.clear();
      return *this;
    }
    
  private:
    
    /// The projector for the full final state.
    FinalState* _fsproj;
    
    /// The final-state particles.
    VetoDetails _vetoCodes;

    
  private:
    
    /// Hide the assignment operator.
    VetoedFinalState& operator=(const VetoedFinalState&);
    
  };

  
}


#endif
