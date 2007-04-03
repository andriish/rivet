// -*- C++ -*-
#ifndef RIVET_DISLepton_H
#define RIVET_DISLepton_H

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/Event/Particle.hh"
#include "Rivet/Tools/Event/Event.hh"

namespace Rivet {

  /// This class projects out the incoming and outgoing leptons in a DIS
  /// event. The incoming lepton is assumed to be along the positive z-axis.
  class DISLepton: public Projection {
    
  public:
    
    /**
     * The default constructor. Must specify the incoming and outgoing
     * PDG codes of the leptons to project.  If \a inid is the
     * anti-particle of \a outid, either a scattered lepton or
     * anti-lepton is searched for. Must also specify a Beam projection
     * object which is assumed to live thoughout the run.
     */
    inline DISLepton(Beam & beamp, long inid, long outid)
      : beams(&beamp), idin(inid), idout(outid) {
      info.declareParameter("BeamA", inid);
      info.declareParameter("DISLepton", outid, "none",
                            "The type of the scattered lepton.");
    }
    
    
  public:
    /// Return the name of the projection
    inline string name() const {
      return "DISLepton";
    }
    
  protected:
    
    /// Perform the projection operation on the supplied event.
    virtual void project(const Event & e);
    
    /// Compare with other projections.
    virtual int compare(const Projection & p) const;
    
  public:
    
    /// The incoming lepton.
    inline const Particle & in() const { return incoming; }
    
    /// The outgoing lepton.
    inline const Particle & out() const { return outgoing; }
    
    /// Return the RivetInfo object of this Projection.
    virtual RivetInfo getInfo() const;
    
  private:
    
    /// The Beam projector object defining the incoming beam particles.
    Beam* beams;
    
    /// The PDG id of the incoming lepton.
    long idin;
    
    /// The PDG id of the outcoming lepton.
    long idout;
    
    /// The incoming lepton.
    Particle incoming;
    
    /// The incoming lepton.
    Particle outgoing;
    
  private:
    
    /// Hide the assignment operator.
    DISLepton & operator=(const DISLepton &);
    
  };
  
}


#endif /* RIVET_DISLepton_H */
