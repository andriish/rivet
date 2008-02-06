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
    Multiplicity(FinalState& fsp)
      : _fsproj(fsp), _totalMult(0), _hadMult(0)
    { 
      addProjection(fsp);
    }

    ~Multiplicity() {
      getLog() << Log::TRACE << "Destroying " << getName() << " at " << this << endl;
    }

  public:
    /// Return the name of the projection
    string getName() const {
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
    const unsigned int totalMultiplicity() const { return _totalMult; }

    /// Hadron multiplicity
    const unsigned int hadronMultiplicity() const { return _hadMult; }
    //@ }

  private:

    /// The FinalState projection used by this projection
    FinalState& _fsproj;

    /// Total multiplicity.
    unsigned int _totalMult;

    /// Hadronic multiplicity.
    unsigned int _hadMult;
  };

}

#endif
