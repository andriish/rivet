// -*- C++ -*-
#ifndef RIVET_Multiplicity_HH
#define RIVET_Multiplicity_HH
// Declaration of the Multiplicity class.

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet {

  /// Project out all final-state particles in an event.
  class Multiplicity : public Projection {

  public:

    /// Constructor. The provided FinalState projection must live throughout the run.
    inline Multiplicity(FinalState& fsp)
      : totalMult_(0), totalChMult_(0), totalUnchMult_(0),
        hadMult_(0), hadChMult_(0), hadUnchMult_(0), fsproj(&fsp)
    { }

  public:
    /// Return the name of the projection
    inline string name() const {
      return "Multiplicity";
    }

  protected:

    /// Perform the projection on the Event.
    void project(const Event & e);

    /// Compare projections.
    int compare(const Projection & p) const;

  public:

    /// @name Access the projected multiplicities.
    //@ {
    /// Total multiplicity
    inline const unsigned int totalMultiplicity() const { return totalMult_; }

    /// Charged multiplicity
    inline const unsigned int totalChargedMultiplicity() const { return totalChMult_; }

    /// Uncharged multiplicity
    inline const unsigned int totalUnchargedMultiplicity() const { return totalUnchMult_; }

    /// Hadron multiplicity
    inline const unsigned int hadronMultiplicity() const { return hadMult_; }

    /// Hadronic charged multiplicity
    inline const unsigned int hadronChargedMultiplicity() const { return hadChMult_; }

    /// Hadronic uncharged multiplicity
    inline const unsigned int hadronUnchargedMultiplicity() const { return hadUnchMult_; }
    //@ }

    
    /// Return the RivetInfo object of this Projection.
    //virtual RivetInfo getInfo() const;

  private:

    /// Total multiplicities
    unsigned int totalMult_, totalChMult_, totalUnchMult_;

    /// Hadronic multiplicities
    unsigned int hadMult_, hadChMult_, hadUnchMult_;

    /// The FinalState projection used by this projection
    FinalState * fsproj;

  private:

    /// Hide the assignment operator.
    Multiplicity & operator=(const Multiplicity &);

  };

}

#endif /* RIVET_Multiplicity_HH */
