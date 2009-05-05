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


  /// Chain together different projections as convenience for finding Z's
  /// from two leptons in the final state
  class ZFinder : public FinalState {

  public:
    
    /// @name Constructors
    //@{

    /// Constructor taking a FinalState and type of the leptons, mass window,
    /// and maximum dR of photons around leptons to take into account for Z
    /// reconstruction.
    ZFinder(const FinalState& fs,
            PdgId pid,
            double m2_min, double m2_max,
            double dRmax);


    /// Constructor taking single eta/pT bounds and type of the leptons, mass
    /// window, and maximum dR of photons around leptons to take into account
    /// for Z reconstruction.
    ZFinder(double etaMin, double etaMax,
            double pTmin,
            PdgId pid,
            double m2_min, double m2_max,
            double dRmax);


    /// Constructor taking multiple eta/pT bounds and type of the leptons, mass
    /// window, and maximum dR of photons around leptons to take into account
    /// for Z reconstruction.
    ZFinder(const std::vector<std::pair<double, double> >& etaRanges,
            double pTmin,
            PdgId pid,
            double m2_min, const double m2_max,
            double dRmax);


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
               double dRmax);

    /// Common implementation of constructor operation, taking FS.
    void _init(const FinalState& fs,
               PdgId pid,
               double m2_min, double m2_max,
               double dRmax);

  };

  
}


#endif
