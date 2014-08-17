// -*- C++ -*-
#ifndef RIVET_FinalState_HH
#define RIVET_FinalState_HH

#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  /// @brief Project out all final-state particles in an event.
  ///
  /// Probably the most important projection in Rivet!
  class FinalState : public ParticleFinder {
  public:

    /// @name Standard constructors and destructors.
    //@{

    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    FinalState(double mineta=-MAXDOUBLE, double maxeta=MAXDOUBLE, double minpt=0.0*GeV);

    /// A constructor which allows to specify multiple eta ranges
    /// and the min \f$ p_T \f$ (in GeV).
    FinalState(const vector<pair<double, double> >& etaRanges, double minpt=0.0*GeV);

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
