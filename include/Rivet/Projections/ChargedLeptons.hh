// -*- C++ -*-
#ifndef RIVET_ChargedLeptons_H
#define RIVET_ChargedLeptons_H

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet {

  /// Project out charged final-state leptons (i.e. electrons and muons, unless you set taus stable!).
  class ChargedLeptons: public Projection {
    
  public:
    
    /// Constructor. The provided FinalState projection must live throughout the run.
    inline ChargedLeptons(FinalState& fsp)
      : fsproj(&fsp)
    { 
      addProjection(fsp);
    }
    
  public:
    /// Return the name of the projection
    inline string getName() const {
      return "ChargedLeptons";
    }
    
  protected:
    
    /// Apply the projection to the event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection & p) const;
    
  public:
    
    /// Access the projected leptons.
    inline const ParticleVector& chargedLeptons() const { return _theChargedLeptons; }
    
  private:
        
    /// The FinalState projection used by this projection
    FinalState * fsproj;

    /// The leptons
    ParticleVector _theChargedLeptons;
    
  private:
    
    /// Hide the assignment operator.
    ChargedLeptons & operator=(const ChargedLeptons &);
    
  };
  
}


#endif
