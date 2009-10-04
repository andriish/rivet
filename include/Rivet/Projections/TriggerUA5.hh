// -*- C++ -*-
#ifndef RIVET_TriggerUA5_HH
#define RIVET_TriggerUA5_HH

#include "Rivet/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// Project out the incoming beams
  class TriggerUA5 : public Projection {
  public:
    
    /// Default constructor.
    TriggerUA5() { 
      setName("TriggerUA5");

      addProjection(Beam(), "Beam");
      addProjection(ChargedFinalState(-3.5, 3.5), "CFS");

      _n_plus = 0;
      _n_minus = 0;
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new TriggerUA5(*this);
    }


  public:

    /// The trigger result for non-single diffractive (2 arm) trigger
    const bool sdDecision() const {
      return _decision_sd;
    }

    /// The trigger result for non-single diffractive (2 arm) trigger
    const bool nsdDecision() const {
      return _decision_nsd;
    }

    /// The trigger result
    const bool samebeams() const {
      return _samebeams;
    }

    /// Number of hits in <-,+> eta hodoscopes
    pair<unsigned int, unsigned int> numHits() {
      return make_pair(_n_plus, _n_minus);
    }

    /// Project on to the Event
    void project(const Event& evt);


  protected:

    /// Compare with other projections.
    virtual int compare(const Projection& p) const {
      return PCmp::EQUIVALENT;
    }


  private:

    /// The min bias trigger decisions
    bool _decision_sd, _decision_nsd;

    /// Is it a pp collision?
    bool _samebeams;

    /// Number of hits in hodoscopes
    unsigned int _n_plus, _n_minus;

  };


}

#endif
