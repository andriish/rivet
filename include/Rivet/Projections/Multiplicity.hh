// -*- C++ -*-
#ifndef RIVET_Multiplicity_HH
#define RIVET_Multiplicity_HH
// Declaration of the Multiplicity class.

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Event/Particle.hh"
#include "Rivet/Tools/Event/Event.hh"


namespace Rivet {

  /**
   * Project out all final-state particles in an event.
   */
  class Multiplicity : public Projection {

  public:

    /// @name Standard constructors and destructors.  @{ / Default
    //constructor. Must specify a FinalState projection which is
    //assumed to live throughout the run.
    inline Multiplicity(FinalState& fsp)
      : totalMult_(0), totalChMult_(0), totalUnchMult_(0),
        hadMult_(0), hadChMult_(0), hadUnchMult_(0), fsproj(&fsp)
    { }

    /// Copy constructor.
    inline Multiplicity(const Multiplicity & x)
      : Projection(x), 
        totalMult_(x.totalMult_), totalChMult_(x.totalChMult_),
        totalUnchMult_(x.totalUnchMult_), hadMult_(x.hadMult_),
        hadChMult_(x.hadChMult_), hadUnchMult_(x.hadUnchMult_),
        fsproj(fsproj)
    { }

    /// Destructor.
    virtual ~Multiplicity() { }
    //@}

  public:
    /// Return the name of the projection
    inline string name() const {
      return "Multiplicity";
    }

  protected:

    /// Perform the projection on the Event: only to be called by 
    /// Event::applyProjection(Projection &).
    void project(const Event & e);

    /// This function defines a unique ordering between different 
    /// Projection objects of the same class. Should only be called from 
    /// operator<(const Projection &).
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

  /**
   * Return the RivetInfo object of this Projection. Derived classes
   * should re-implement this function to return the combined
   * RivetInfo object of this and of any other Projection upon which
   * this depends.
   */
  virtual RivetInfo getInfo() const;

  private:

    /// Total multiplicities
    unsigned int totalMult_, totalChMult_, totalUnchMult_;

    /// Hadronic multiplicities
    unsigned int hadMult_, hadChMult_, hadUnchMult_;

    /// The FinalState projection used by this projection
    FinalState * fsproj;

  private:

    /**
     * The assignment operator is private and must never be called.
     * In fact, it shouldn't even be implemented.
     */
    Multiplicity & operator=(const Multiplicity &);

  };

}

#endif /* RIVET_Multiplicity_HH */
