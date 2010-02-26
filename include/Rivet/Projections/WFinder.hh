// -*- C++ -*-
#ifndef RIVET_WFinder_HH
#define RIVET_WFinder_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Chain together different projections as convenience for finding W's
  /// from two leptons in the final state
  class WFinder : public FinalState {
  public:
 
    /// @name Constructors
    //@{

    /// Constructor taking a FinalState and type of the charged lepton, mass window,
    /// and maximum dR of photons around the charged lepton to take into account for W
    /// reconstruction.
    WFinder(const ChargedFinalState& fs_l,
            PdgId pid,
            double m2_min, double m2_max,
            double missingET,
            double dRmax);


    /// Constructor taking single eta/pT bounds and type of the charged lepton, mass
    /// window, and maximum dR of photons around the charged lepton to take into account
    /// for W reconstruction.
    WFinder(double etaMin, double etaMax,
            double pTmin,
            PdgId pid,
            double m2_min, double m2_max,
            double missingET,
            double dRmax);


    /// Constructor taking multiple eta/pT bounds and type of the charged lepton, mass
    /// window, and maximum dR of photons around the charged lepton to take into account
    /// for W reconstruction.
    WFinder(const std::vector<std::pair<double, double> >& etaRanges,
            double pTmin,
            PdgId pid,
            double m2_min, const double m2_max,
            double missingET,
            double dRmax);


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new WFinder(*this);
    }
    //@}


    /// Access to the remaining particles, after the W and clustered photons
    /// have been removed from the full final state
    /// (e.g. for running a jet finder on it)
    const FinalState& remainingFinalState() const;

    /// Access to the W constituent leptons and photons
    const FinalState& constituentsFinalState() const;

    /// Access to the W constituent leptons
    const FinalState& constituentLeptonsFinalState() const;


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
               double missingET,
               double dRmax);

    /// Common implementation of constructor operation, taking FS.
    void _init(const ChargedFinalState& fs_l,
               PdgId pid,
               double m2_min, double m2_max,
               double missingET,
               double dRmax);


  private:

    // Mass range
    double _m2_min, _m2_max;

    // Missing ET cut
    double _etMiss;

  };


}


#endif
