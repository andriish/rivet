// -*- C++ -*-
#ifndef RIVET_LeptonClusters_HH
#define RIVET_LeptonClusters_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  class ClusteredLepton : public Particle {
  public:

    ClusteredLepton(Particle lepton) :
      Particle(lepton.pdgId(), lepton.momentum()),
      _constituentLepton(lepton) {}

    void addPhoton(const Particle& p, bool cluster) {
      _constituentPhotons.push_back(p);
      if (cluster) setMomentum(momentum() + p.momentum());
    }

    const Particle& constituentLepton() const { return _constituentLepton; }
    const ParticleVector& constituentPhotons() const { return _constituentPhotons; }

  private:

    ParticleVector _constituentPhotons;
    Particle _constituentLepton;
  };


  /// @brief Cluster photons from a given FS to all charged particles (typically
  /// leptons) from signal and store the original charged particles and photons
  /// as particles() while the newly created clustered lepton objects are
  /// accessible as clusteredLeptons()
  /// @brief Use LeptonClustersConstituents projection to cluster all photons to
  /// leptons. This projection here makes the clustered objects available
  /// in a FinalState manner, i.e. as particles(). The given pT and eta cuts are
  /// applied to those objects. The underlying (original) leptons and photons
  /// are available by means of the constituentsFinalState() method.
  class LeptonClusters : public FinalState {

  public:

    LeptonClusters(const FinalState& photons, const FinalState& signal,
                   double dRmax, bool cluster,
                   const std::vector<std::pair<double, double> >& etaRanges,
                   double pTmin);

    virtual const Projection* clone() const {
      return new LeptonClusters(*this);
    }

    const vector<ClusteredLepton>& clusteredLeptons() const { return _clusteredLeptons; }

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

    /// Container which stores the clustered lepton objects
    vector<ClusteredLepton> _clusteredLeptons;
  };




}


#endif
