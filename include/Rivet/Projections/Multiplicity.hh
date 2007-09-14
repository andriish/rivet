// -*- C++ -*-
#ifndef RIVET_Multiplicity_HH
#define RIVET_Multiplicity_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet {

  /// Count the final-state particles in an event.
  class Multiplicity : public Projection {

  public:

    /// Constructor. The provided FinalState projection must live throughout the run.
    inline Multiplicity(FinalState& fsp)
      : _totalMult(0), _totalChMult(0), _totalUnchMult(0),
        _hadMult(0), _hadChMult(0), _hadUnchMult(0), 
        _fsproj(&fsp)
    { 
      addProjection(fsp);
    }

  public:
    /// Return the name of the projection
    inline string getName() const {
      return "Multiplicity";
    }

  protected:

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;

  public:

    /// @name Access the projected multiplicities.
    //@ {
    /// Total multiplicity
    inline const unsigned int totalMultiplicity() const { return _totalMult; }

    /// Charged multiplicity
    inline const unsigned int totalChargedMultiplicity() const { return _totalChMult; }

    /// Uncharged multiplicity
    inline const unsigned int totalUnchargedMultiplicity() const { return _totalUnchMult; }

    /// Hadron multiplicity
    inline const unsigned int hadronMultiplicity() const { return _hadMult; }

    /// Hadronic charged multiplicity
    inline const unsigned int hadronChargedMultiplicity() const { return _hadChMult; }

    /// Hadronic uncharged multiplicity
    inline const unsigned int hadronUnchargedMultiplicity() const { return _hadUnchMult; }
    //@ }

  private:

    /// Total multiplicities
    unsigned int _totalMult, _totalChMult, _totalUnchMult;

    /// Hadronic multiplicities
    unsigned int _hadMult, _hadChMult, _hadUnchMult;

    /// The FinalState projection used by this projection
    FinalState* _fsproj;

  };

}

#endif
