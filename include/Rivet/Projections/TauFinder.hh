// -*- C++ -*-
#ifndef RIVET_TauFinder_HH
#define RIVET_TauFinder_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Convenience finder of unstable taus
  /// @todo Inherit directly from ParticleFinder, not FinalState
  class TauFinder : public FinalState {
  public:

    enum DecayType { ANY=0, LEPTONIC=1, HADRONIC };

    /// @name Constructors
    //@{

    /// @todo Why accept a FinalState? Find taus which decay to particles in this FS? Document the logic.
    TauFinder(const FinalState& inputfs, DecayType decaytype=ANY) {
      _init(inputfs, decaytype);
    }

    TauFinder(DecayType decaytype=ANY) {
      _init(UnstableFinalState(), decaytype);
    }


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new TauFinder(*this);
    }

    //@}


    /// Access to the found bosons
    ///
    /// @note Currently either 0 or 1 boson can be found.
    const Particles& taus() const { return _taus; }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    //int compare(const Projection& p) const;


  public:

    /// Clear the projection
    void clear() {
      _theParticles.clear();
      _taus.clear();
    }



  private:

    /// Common implementation of constructor operation, taking FS params.
    void _init(const FinalState& inputfs, DecayType decaytype);

    /// List of found taus
    /// @todo Fill _theParticles instead, when inheriting from ParticleFinder?
    Particles _taus;

  };


}


#endif
