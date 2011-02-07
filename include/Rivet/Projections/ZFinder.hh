// -*- C++ -*-
#ifndef RIVET_ZFinder_HH
#define RIVET_ZFinder_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Convenience finder of leptonically decaying Zs
  ///
  /// Chain together different projections as convenience for finding Z's
  /// from two leptons in the final state, including photon clustering.
  class ZFinder : public FinalState {

  public:

    /// @name Constructors
    //@{

    /// Constructor taking a FinalState and type of the leptons, mass window,
    /// and maximum dR of photons around leptons to take into account for Z
    /// reconstruction (both for clustering to leptons, and exclusion from jets).
    ZFinder(const FinalState& fs,
            PdgId pid,
            double m2_min, double m2_max,
            double dRmax_clustering,
            double dRmax_exclusion);


    /// Constructor taking single eta/pT bounds and type of the leptons, mass
    /// window, and maximum dR of photons around leptons to take into account
    /// for Z reconstruction (both for clustering to leptons, and exclusion from jets).
    ZFinder(double etaMin, double etaMax,
            double pTmin,
            PdgId pid,
            double m2_min, double m2_max,
            double dRmax_clustering,
            double dRmax_exclusion);


    /// Constructor taking multiple eta/pT bounds and type of the leptons, mass
    /// window, and maximum dR of photons around leptons to take into account
    /// for Z reconstruction (both for clustering to leptons, and exclusion from jets).
    ZFinder(const std::vector<std::pair<double, double> >& etaRanges,
            double pTmin,
            PdgId pid,
            double m2_min, const double m2_max,
            double dRmax_clustering,
            double dRmax_exclusion);


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new ZFinder(*this);
    }
    //@}


    /// Access to the remaining particles, after the Z and clustered photons
    /// have been removed from the full final state
    /// (e.g. for running a jet finder on it)
    const FinalState& remainingFinalState() const;

    /// Access to the Z constituent leptons final state
    /// (e.g. for more fine-grained cuts on the leptons)
    const FinalState& constituentsFinalState() const;

    /// Access to the photons which have been clustered to the leptons
    const FinalState& clusteredPhotonsFinalState() const;

  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  private:
    /// Common implementation of constructor operation, taking FS params.
    void _init(const std::vector<std::pair<double, double> >& etaRanges,
               double pTmin,  PdgId pid,
               double m2_min, double m2_max,
               double dRmax_clustering,
               double dRmax_exclusion);

    /// Common implementation of constructor operation, taking FS.
    void _init(const FinalState& fs,
               PdgId pid,
               double m2_min, double m2_max,
               double dRmax_clustering,
               double dRmax_exclusion);

  };


}


#endif
