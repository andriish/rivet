// -*- C++ -*-
#ifndef RIVET_TriggerCDFRun2_HH
#define RIVET_TriggerCDFRun2_HH

#include "Rivet/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Access to the min bias triggers used by CDF in Run 0 and Run 1
  class TriggerCDFRun2 : public Projection {
  public:

    /// Default constructor.
    TriggerCDFRun2() {
      setName("TriggerCDFRun2");

      declare(ChargedFinalState(Cuts::etaIn(-4.7, 4.7)), "CFS");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(TriggerCDFRun2);


    /// The trigger result
    bool minBiasDecision() const {
      return _decision_mb;
    }

    /// Project on to the Event
    void project(const Event& evt);


  protected:

    /// Compare with other projections.
    virtual CmpState compare(const Projection&) const {
      return CmpState::EQ;
    }


  private:

    /// The min bias trigger decision
    bool _decision_mb;

  };


}

#endif
