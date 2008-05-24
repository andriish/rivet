// -*- C++ -*-
#ifndef RIVET_DISLepton_HH
#define RIVET_DISLepton_HH

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// This class projects out the incoming and outgoing leptons in a DIS
  /// event. The incoming lepton is assumed to be along the positive z-axis.
  class DISLepton : public Projection {
    
  public:
    
    /// @name Constructors.
    //@{
    /// Specify the incoming and outgoing PDG codes of the leptons to project.
    /// If \a inid is an anti-particle and \a outid a particle, or vice versa,
    /// either a scattered lepton or anti-lepton is searched for. This version
    /// also specifies a FinalState projection.
    DISLepton(const FinalState& fsp, const ParticleName& inid, const ParticleName& outid)
      : _idin(inid), _idout(outid)
    {
      setName("DISLepton");
      addBeamPair(inid, ANY);
      addProjection(Beam(), "Beam");
      addProjection(fsp, "FS");
    }
    
    /// Specify the incoming and outgoing PDG codes of the leptons to project.
    /// If \a inid is an anti-particle and \a outid a particle, or vice versa,
    /// either a scattered lepton or anti-lepton is searched for. This version
    /// uses a default-constructed basic FinalState projection.
    DISLepton(const ParticleName& inid, const ParticleName& outid)
      : _idin(inid), _idout(outid)
    {
      setName("DISLepton");
      addBeamPair(inid, ANY);
      addProjection(Beam(), "Beam");
      addProjection(FinalState(), "FS");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new DISLepton(*this);
    }
    //@}
    
    
  protected:
    
    /// Perform the projection operation on the supplied event.
    virtual void project(const Event& e);
    
    /// Compare with other projections.
    virtual int compare(const Projection& p) const;
    
  public:
    
    /// The incoming lepton.
    const Particle& in() const { return _incoming; }
    
    /// The outgoing lepton.
    const Particle& out() const { return _outgoing; }
    
  private:
    
    /// The PDG id of the incoming lepton.
    long _idin;
    
    /// The PDG id of the outcoming lepton.
    long _idout;
    
    /// The incoming lepton.
    Particle _incoming;
    
    /// The incoming lepton.
    Particle _outgoing;
        
  };
  
}


#endif
