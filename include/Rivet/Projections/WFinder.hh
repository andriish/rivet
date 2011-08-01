// -*- C++ -*-
#ifndef RIVET_WFinder_HH
#define RIVET_WFinder_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/LeptonClusters.hh"

namespace Rivet {


  /// @brief Convenience finder of leptonically decaying Ws
  ///
  /// Chain together different projections as convenience for finding W's
  /// from two leptons in the final state, including photon clustering.
  class WFinder : public FinalState {
  public:

    /// @name Constructors
    //@{

    // /// Constructor taking a FinalState and type of the charged lepton, mass window,
    // /// and maximum dR of photons around the charged lepton to take into account for W
    // /// reconstruction.
    // WFinder(const ChargedFinalState& fs_l,
    //         PdgId pid,
    //         double minmass, double maxmass,
    //         double missingET,
    //         double dRmax, bool clusterPhotons, bool excludePhotonsFromRFS);


    /// Constructor taking single eta/pT bounds and type of the charged lepton, transverse mass
    /// window, and maximum dR of photons around the charged lepton to take into account
    /// for W reconstruction.
    WFinder(double etaMin, double etaMax,
            double pTmin,
            PdgId pid,
            double minmass, double maxmass,
            double missingET,
            double dRmax, bool clusterPhotons=true, bool trackPhotons=false,
            double masstarget=80.4,
            bool useTransverseMass=false);


    /// Constructor taking multiple eta/pT bounds and type of the charged lepton, transverse mass
    /// window, and maximum dR of photons around the charged lepton to take into account
    /// for W reconstruction.
    WFinder(const std::vector<std::pair<double, double> >& etaRanges,
            double pTmin,
            PdgId pid,
            double minmass, const double maxmass,
            double missingET,
            double dRmax, bool clusterPhotons=true, bool trackPhotons=false,
            double masstarget=80.4,
            bool useTransverseMass=false);


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new WFinder(*this);
    }
    //@}


    /// Access to the found bosons (currently either 0 or 1)
    const ParticleVector& bosons() const { return _bosons; }

    /// Access to the W constituent clustered leptons (currently either of
    /// size 0 if no boson was found or 1 if one boson was found)
    const vector<Particle>& constituentLeptons() const { return _constituentLeptons; }

    /// Access to the W constituent neutrinos (currently either of size 0 if no
    /// boson was found or 1 if one boson was found)
    const vector<Particle>& constituentNeutrinos() const { return _constituentNeutrinos; }

    /// Access to the remaining particles, after the W and clustered photons
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
      _constituentLeptons.clear();
      _constituentNeutrinos.clear();
    }


  private:

    /// Common implementation of constructor operation, taking FS params.
    void _init(const std::vector<std::pair<double, double> >& etaRanges,
               double pTmin,  PdgId pid,
               double minmass, double maxmass,
               double missingET,
               double dRmax, bool clusterPhotons, bool trackPhotons,
               double masstarget,
               bool useTransverseMass);


  private:

    /// Transverse mass cuts
    double _minmass, _maxmass, _masstarget;
    bool _useTransverseMass;

    /// Missing ET cut
    double _etMiss;

    /// Switch for tracking of photons (whether to add them to _theParticles)
    /// This is relevant when the ZFinder::_theParticles are to be excluded
    /// from e.g. the input to a jet finder, to specify whether the clustered
    /// photons are to be excluded as well.
    /// (Yes, some experiments make a difference between clusterPhotons and
    /// trackPhotons!)
    bool _trackPhotons;

    /// Lepton flavour
    PdgId _pid;

    /// Neutrino flavour
    PdgId _nu_pid;

    /// list of found bosons (currently either 0 or 1)
    ParticleVector _bosons;

    /// Constituent leptons (currently either 0 or 1)
    ParticleVector _constituentLeptons;

    /// Constituent neutrinos (currently either 0 or 1)
    ParticleVector _constituentNeutrinos;

  };


}


#endif
