// -*- C++ -*-
#ifndef RIVET_TRIGGERPROJECTION_HH
#define RIVET_TRIGGERPROJECTION_HH

#include "Rivet/Projection.hh"

namespace Rivet {

  /** @brief Base class for projections returning a bool corresponding
      to a trigger.

      @author Leif Lönnblad

      Project an event down to a single true or false value accessible
      through the operator() function, where true means that the event
      has passed some trigger criterion.

  */
  class TriggerProjection: public Projection {

  public:

    /// The default constructor.
    TriggerProjection() : _passed(true) {
      setName("TriggerProjection");
    }
    virtual ~TriggerProjection() {}

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(TriggerProjection);

    /// Return true if the event has passed some trigger or selection
    /// criteria.
    bool operator()() const {
      return _passed;
    }

  protected:

    virtual void project(const Event& e) {
      pass();
    }

    /// Indicate that the event has passed the trigger.
    void pass() {
      _passed = true;
    }

    /// Compare projections
    virtual CmpState compare(const Projection&) const {
      return CmpState::EQ;
    }

    /// Indicate that the event has failed the trigger.
    void fail() {
      _passed = false;
    }

  protected:

    bool _passed;

  };

}

#endif
