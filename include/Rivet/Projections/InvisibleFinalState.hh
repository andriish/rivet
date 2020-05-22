// -*- C++ -*-
#ifndef RIVET_InvisibleFinalState_HH
#define RIVET_InvisibleFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Final state modifier excluding particles which are experimentally visible
  class InvisibleFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    /// Constructor with specific FinalState.
    InvisibleFinalState(bool requirepromptness=false)
      : _requirePromptness(requirepromptness)
    {
      setName("InvisibleFinalState");
      declare(FinalState(), "FS");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(InvisibleFinalState);

    //@}

    /// Require accepted particles to be prompt
    void requirePromptness(bool acc=true) { _requirePromptness = acc; }

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;

    private:

      bool _requirePromptness;
  };


}

#endif
