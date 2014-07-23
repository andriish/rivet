// -*- C++ -*-
#ifndef RIVET_TauFinder_HH
#define RIVET_TauFinder_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Convenience finder of leptonically decaying Ws
  ///
  /// Chain together different projections as convenience for finding W's
  /// from two leptons in the final state, including photon clustering.
  class TauFinder : public FinalState {
  public:

    enum DecayType { ANY=0, LEPTONIC=1, HADRONIC };

    /// @name Constructors
    //@{

    ///
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



    /// Useful for e.g. input to a jet finder
    //const FinalState& remainingFinalState() const;



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

    /// list of found bosons (currently either 0 or 1)
    Particles _taus;

  };


}


#endif

