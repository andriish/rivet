// -*- C++ -*-
#ifndef RIVET_HeavyHadrons_HH
#define RIVET_HeavyHadrons_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Project out the last pre-decay b and c hadrons.
  ///
  /// This currently defines a c-hadron as one which contains a @a c quark and
  /// @a{not} a @a b quark.
  ///
  /// @todo This assumes that the heavy hadrons are unstable... should we also look for stable ones in case the decays are disabled?
  class HeavyHadrons : public FinalState {
  public:

    /// @name Constructors and destructors.
    //@{

    /// Constructor with specification of the minimum and maximum pseudorapidity
    /// \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    HeavyHadrons(double mineta = -MAXRAPIDITY,
                 double maxeta =  MAXRAPIDITY,
                 double minpt  =  0.0*GeV) {
      setName("HeavyHadrons");
      addProjection(UnstableFinalState(mineta, maxeta, minpt), "UFS");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new HeavyHadrons(*this);
    }

    //@}


    /// @name Particle accessors
    //@{

    /// Get the pre-decay b hadrons
    const Particles& bHadrons() const {
      return _theBs;
    }

    /// Get the pre-decay c hadrons
    const Particles& cHadrons() const {
      return _theCs;
    }

    //@}


  protected:

    /// Apply the projection to the event.
    virtual void project(const Event& e);

    /// b and c hadron containers
    Particles _theBs, _theCs;

  };


}


#endif
