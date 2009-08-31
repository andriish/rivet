// -*- C++ -*-
#ifndef RIVET_LeadingParticlesFinalState_HH
#define RIVET_LeadingParticlesFinalState_HH

#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

  class Particle;

  /// Project only the highest-pT occurrences of FS particles with the specified PDG IDs.
  class LeadingParticlesFinalState : public FinalState {

  public:

    /// Constructor: the supplied FinalState projection is assumed to live through the run.
    /// @todo Why specify the rap & pT cuts again?
    LeadingParticlesFinalState(FinalState & fsp, double mineta = -MAXRAPIDITY, double maxeta = MAXRAPIDITY, double minpt = 0.0 * GeV)
  :  FinalState(mineta, maxeta, minpt) {
      setName("LeadingParticlesFinalState");
      addProjection(fsp, "FS");
    }

    /// Clone on the heap. 
    virtual const Projection *clone() const {
      return new LeadingParticlesFinalState(*this);
    }

    /// Add a particle ID to the list of leading particles selected 
    LeadingParticlesFinalState& addParticleId(long id) {
      _ids.insert(id);
      return *this;
    } 

    /// Add a particle ID to the list of leading particles selected 
    LeadingParticlesFinalState& addParticleIdPair(long id) {
      _ids.insert(id);
      _ids.insert(-id);
      return *this;
    } 
    
  protected:

    /// Apply the projection on the supplied event.
    void project(const Event & e);

    /// Compare projections.
    int compare(const Projection & p) const;

    /// check if the particle's id is in the list 
    bool inList(const Particle & particle) const;

  private:

    /// IDs of the leading particles to be selected
    std::set < long >_ids;

  };

}

#endif
