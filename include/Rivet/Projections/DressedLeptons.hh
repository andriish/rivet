// -*- C++ -*-
#ifndef RIVET_DressedLeptons_HH
#define RIVET_DressedLeptons_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Cuts.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// A charged lepton meta-particle created by clustering photons close to the bare lepton
  class DressedLepton : public Particle {
  public:

    DressedLepton(const Particle& lepton) :
      Particle(lepton.pdgId(), lepton.momentum()),
      _constituentLepton(lepton) {}

    void addPhoton(const Particle& p, bool cluster) {
      _constituentPhotons.push_back(p);
      if (cluster) setMomentum(momentum() + p.momentum());
    }

    const Particle& constituentLepton() const { return _constituentLepton; }
    const Particles& constituentPhotons() const { return _constituentPhotons; }

  private:

    Particles _constituentPhotons;
    Particle _constituentLepton;
  };


  /// @brief Cluster photons from a given FS to all charged particles (typically leptons)
  ///
  /// This stores the original charged particles and photons as particles()
  /// while the newly created clustered lepton objects are accessible as
  /// clusteredLeptons().
  class DressedLeptons : public FinalState {
  public:

    /// Constructor with a single eta range
    DressedLeptons(const FinalState& photons, const FinalState& signal,
                   double dRmax, bool cluster,
                   double etaMin, double etaMax,
                   double pTmin, bool useDecayPhotons=false);

    /// Constructor with multiple eta ranges
    DressedLeptons(const FinalState& photons, const FinalState& signal,
                   double dRmax, bool cluster, Cut c,
                   bool useDecayPhotons=false);


    /// Clone this projection
    virtual const Projection* clone() const {
      return new DressedLeptons(*this);
    }

    /// Retrieve the dressed leptons
    const vector<DressedLepton>& dressedLeptons() const { return _clusteredLeptons; }

    /// Retrieve the dressed leptons (synonym)
    /// @deprecated Use dressedLeptons()
    const vector<DressedLepton>& clusteredLeptons() const { return _clusteredLeptons; }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  private:

    /// Maximum cone radius to find photons in
    double _dRmax;
    /// Whether to actually add the photon momenta to clusteredLeptons
    bool _cluster;
    /// Whether to include photons from hadron (particularly pi0) decays
    bool _fromDecay;

    /// Container which stores the clustered lepton objects
    vector<DressedLepton> _clusteredLeptons;

  };



}


#endif
