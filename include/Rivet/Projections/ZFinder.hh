// -*- C++ -*-
#ifndef RIVET_ZFinder_HH
#define RIVET_ZFinder_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeptonClusters.hh"

namespace Rivet {


  /// @brief Convenience finder of leptonically decaying Zs
  ///
  /// Chain together different projections as convenience for finding Z's
  /// from two leptons in the final state, including photon clustering.
  class ZFinder : public FinalState {

  public:

    /// @name Constructors
    //@{

    /// Constructor taking single eta/pT bounds
    /// @param inputfs Input final state
    /// @param etaMin,etaMax,pTmin lepton cuts
    /// @param pid type of the leptons
    /// @param minmass,maxmass mass window
    /// @param dRmax maximum dR of photons around leptons to take into account
    ///  for Z reconstruction (only relevant if one of the following are true)
    /// @param clusterPhotons whether such photons are supposed to be
    ///  clustered to the lepton objects and thus Z mom
    /// @param trackPhotons whether such photons should be added to _theParticles
    ///  (cf. _trackPhotons)
    ZFinder(const FinalState& inputfs,
            Cut cuts,
            PdgId pid,
            double minmass, double maxmass,
            double dRmax, bool clusterPhotons, bool trackPhotons,
            double masstarget=91.2*GeV);

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new ZFinder(*this);
    }
    //@}


    /// Access to the found bosons (currently either 0 or 1)
    const Particles& bosons() const { return _bosons; }

    /// Access to the Z constituent clustered leptons
    /// (e.g. for more fine-grained cuts on the clustered leptons)
    /// The order is going to be: positive charge constituent 1st, negative 2nd
    const vector<Particle>& constituents() const { return _constituents; }

    /// Access to the remaining particles, after the Z and clustered photons
    /// have been removed from the full final state
    /// (e.g. for running a jet finder on it)
    const FinalState& remainingFinalState() const;


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  public:

    /// Clear the projection
    void clear() {
      _theParticles.clear();
      _bosons.clear();
      _constituents.clear();
    }


  private:
    /// Mass cuts to apply to clustered leptons (cf. InvMassFinalState)
    double _minmass, _maxmass, _masstarget;

    /// Switch for tracking of photons (whether to add them to _theParticles)
    /// This is relevant when the ZFinder::_theParticles are to be excluded
    /// from e.g. the input to a jet finder, to specify whether the clustered
    /// photons are to be excluded as well.
    /// (Yes, some experiments make a difference between clusterPhotons and
    /// trackPhotons!)
    bool _trackPhotons;

    /// Lepton flavour
    PdgId _pid;

    /// list of found bosons (currently either 0 or 1)
    Particles _bosons;

    /// Clustered leptons
    vector<Particle> _constituents;

  };


}



#endif
