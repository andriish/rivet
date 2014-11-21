// -*- C++ -*-
#ifndef RIVET_FinalState_HH
#define RIVET_FinalState_HH

#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Cuts.hh"

namespace Rivet {


  /// @brief Project out all final-state particles in an event.
  /// Probably the most important projection in Rivet!
  class FinalState : public ParticleFinder {
  public:

    /// @name Standard constructors and destructors.
    //@{
    /// @deprecated Keep for backwards compatibility for now
    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    FinalState(double mineta,
               double maxeta,
               double minpt = 0.0);

    /// Construction using Cuts object
    FinalState(const Cut & c = Cuts::open());

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new FinalState(*this);
    }

    //@}


    /// Apply the projection to the event.
    virtual void project(const Event& e);

    /// Compare projections.
    virtual int compare(const Projection& p) const;

    /// Decide if a particle is to be accepted or not.
    /// @todo Rename to _accept or acceptFinal?
    virtual bool accept(const Particle& p) const;

  };


}

#endif
