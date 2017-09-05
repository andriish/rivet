// -*- C++ -*-
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  ZFinder::ZFinder(const FinalState& inputfs,
		   const Cut & fsCut,
                   PdgId pid,
                   double minmass, double maxmass,
                   double dRmax,
                   ClusterPhotons clusterPhotons,
                   PhotonTracking trackPhotons,
                   double masstarget)
  {
    setName("ZFinder");

    _minmass = minmass;
    _maxmass = maxmass;
    _masstarget = masstarget;
    _pid = abs(pid);
    _trackPhotons = trackPhotons;

    IdentifiedFinalState bareleptons(inputfs);
    bareleptons.acceptIdPair(_pid);
    const bool doClustering = (clusterPhotons != NOCLUSTER);
    const bool useDecayPhotons = (clusterPhotons == CLUSTERALL);
    DressedLeptons leptons(inputfs, bareleptons, (doClustering ? dRmax : -1.0), fsCut, useDecayPhotons);
    addProjection(leptons, "DressedLeptons");

    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(*this);
    addProjection(remainingFS, "RFS");
  }



  /////////////////////////////////////////////////////



  Particles ZFinder::constituentLeptons() const {
    if (empty()) return Particles();
    return boson().constituents();
    // return boson().constituents(isChargedLepton);
  }


  Particles ZFinder::constituentLeptons(const ParticleSorter& cmp) const {
    return sortBy(constituentLeptons(), cmp);
  }


  const VetoedFinalState& ZFinder::remainingFinalState() const {
    return getProjection<VetoedFinalState>("RFS");
  }


  int ZFinder::compare(const Projection& p) const {
    PCmp LCcmp = mkNamedPCmp(p, "DressedLeptons");
    if (LCcmp != EQUIVALENT) return LCcmp;

    const ZFinder& other = dynamic_cast<const ZFinder&>(p);
    return (cmp(_minmass, other._minmass) ||
            cmp(_maxmass, other._maxmass) ||
            cmp(_pid, other._pid) ||
            cmp(_trackPhotons, other._trackPhotons));
  }


  void ZFinder::project(const Event& e) {
    clear();

    // Get leptons and find an acceptable invariant mass OSSF pair
    const DressedLeptons& leptons = applyProjection<DressedLeptons>(e, "DressedLeptons");
    InvMassFinalState imfs(std::make_pair(_pid, -_pid), _minmass, _maxmass, _masstarget);
    auto dressed = leptons.dressedLeptons();
    Particles tmp(dressed.begin(), dressed.end());
    imfs.calc(tmp);
    if (imfs.particlePairs().size() < 1) return;

    // Assemble a pseudo-Z particle
    const ParticlePair& Zconstituents = imfs.particlePairs().front();
    const Particle& p1(Zconstituents.first), p2(Zconstituents.second);
    const FourMomentum pZ = p1.momentum() + p2.momentum();
    assert(p1.charge3() + p2.charge3() == 0);
    Particle z(PID::Z0BOSON, pZ);

    // Debug printout
    stringstream msg;
    msg << "Z " << pZ << " reconstructed from: \n"
        << "   " << p1.momentum() << " " << p1.pid() << "\n"
        << " + " << p2.momentum() << " " << p2.pid();
    MSG_DEBUG(msg.str());

    // Add (dressed) lepton constituents to the W (skipping photons if requested)
    // Keep the DressedLeptons found by the ZFinder
    const DressedLepton l1 = p1.charge() > 0 ? p1 : p2;
    z.addConstituent(_trackPhotons ? l1 : l1.bareLepton());
    const DressedLepton l2 = p2.charge() < 0 ? p2 : p1;
    z.addConstituent(_trackPhotons ? l2 : l2.bareLepton());

    // Register the completed Z
    _theParticles.push_back(z);
  }


}
