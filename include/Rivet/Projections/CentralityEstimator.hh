// -*- C++ -*-
#ifndef RIVET_CentralityEstimator_HH
#define RIVET_CentralityEstimator_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Event.hh"
#include "HepMC/GenEvent.h"

namespace Rivet {


  /**
    @brief Base class for projections giving the value of an
    observable sensitive to the centrality of a collision.

    @author Leif Lönnblad

    The centrality of a collision is not really an observable, but the
    concept is anyway often used in the heavy ion community as if it
    were just that.

    This base class can be used to provide a an estimator for the
    centrality by projecting down to a single number which then canbe
    used by a CentralityGroup object to select a histogram to be
    filled with another observable depending on centrality percentile.

    The eztimate() should be a non-negative number with large values
    indicating a higher overlap than small ones. A negative value
    indicates that the centrality estimate could not be calculated.

    In the best of all worlds the centrality estimator should be a
    proper hadron-level observable corrected for detector effects,
    however, this base class only returns the inverse of the
    impact_parameter member of the GenHeavyIon object in an GenEvent
    if present and zero otherwise.
   */
  class CentralityEstimator : public Projection {
  public:

    /// Constructor.
    CentralityEstimator()
      : _estimate(-1.0) {}

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(CentralityEstimator);

  protected:

    /// Perform the projection on the Event
    void project(const Event& e) {
      _estimate = -1.0;
      const HepMC::HeavyIon * hi = e.genEvent()->heavy_ion();
      if ( hi ) _estimate = hi->impact_parameter() > 0.0?
                  1.0/hi->impact_parameter(): numeric_limits<double>::max();
    }

    /// Compare projections
    int compare(const Projection& p) const {
      return mkNamedPCmp(p, "CentEst");
    }


  public:

    /// The value of the centrality estimate.
    double estimate() const { return _estimate; }


  protected:

    /// The value of the centrality estimate.
    double _estimate;

  };

}

#endif

