// -*- C++ -*-
#ifndef RIVET_DISLepton_H
#define RIVET_DISLepton_H

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {

  /// This class projects out the incoming and outgoing leptons in a DIS
  /// event. The incoming lepton is assumed to be along the positive z-axis.
  class DISLepton: public Projection {
    
  public:
    
    /// The default constructor. Must specify the incoming and outgoing
    /// PDG codes of the leptons to project.  If \a inid is the
    /// anti-particle of \a outid, either a scattered lepton or
    /// anti-lepton is searched for. Must also specify a Beam projection
    /// object which is assumed to live thoughout the run.
    inline DISLepton(Beam& beamproj, const ParticleName& inid, const ParticleName& outid)
      : _beams(&beamproj), _idin(inid), _idout(outid) {
      _beamPairs.insert(BeamPair(inid, ANY));
      addProjection(beamproj);
    }
    
    
  public:
    /// Return the name of the projection
    inline string getName() const {
      return "DISLepton";
    }
    
  protected:
    
    /// Perform the projection operation on the supplied event.
    virtual void project(const Event & e);
    
    /// Compare with other projections.
    virtual int compare(const Projection & p) const;
    
  public:
    
    /// The incoming lepton.
    inline const Particle & in() const { return _incoming; }
    
    /// The outgoing lepton.
    inline const Particle & out() const { return _outgoing; }
    
  private:
    
    /// The Beam projector object defining the incoming beam particles.
    Beam* _beams;
    
    /// The PDG id of the incoming lepton.
    long _idin;
    
    /// The PDG id of the outcoming lepton.
    long _idout;
    
    /// The incoming lepton.
    Particle _incoming;
    
    /// The incoming lepton.
    Particle _outgoing;
    
  private:
    
    /// Hide the assignment operator.
    DISLepton & operator=(const DISLepton &);
    
  };
  
}


#endif /* RIVET_DISLepton_H */
