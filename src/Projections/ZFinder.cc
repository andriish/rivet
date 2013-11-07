// -*- C++ -*-
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/LeptonClusters.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  ZFinder::ZFinder(const FinalState& inputfs,
		   Cut cuts,
                   PdgId pid,
                   double minmass, double maxmass,
                   double dRmax, bool clusterPhotons, bool trackPhotons,
                   double masstarget) {
    _init(inputfs, cuts, 
	  pid, minmass, maxmass, dRmax, clusterPhotons, trackPhotons, masstarget);
  }

  ZFinder::ZFinder(const FinalState& inputfs,
                   double etaMin, double etaMax,
                   double pTmin,
                   PdgId pid,
                   double minmass, double maxmass,
                   double dRmax, bool clusterPhotons, bool trackPhotons,
                   double masstarget) {
    Cut eta = Range( Cuts::eta, etaMin, etaMax );
    Cut pt  = Cuts::pt >= pTmin;
    _init(inputfs, eta & pt, pid, minmass, maxmass, dRmax, clusterPhotons, trackPhotons, masstarget);
  }


  // ZFinder::ZFinder(const FinalState& inputfs,
  //                  const std::vector<std::pair<double, double> >& etaRanges,
  //                  double pTmin,
  //                  PdgId pid,
  //                  double minmass, const double maxmass,
  //                  double dRmax, bool clusterPhotons, bool trackPhotons,
  //                  double masstarget) {
  //   _init(inputfs, etaRanges, pTmin, pid, minmass, maxmass, dRmax, clusterPhotons, trackPhotons, masstarget);
  // }


  // ZFinder::ZFinder(double etaMin, double etaMax,
  //                  double pTmin,
  //                  PdgId pid,
  //                  double minmass, double maxmass,
  //                  double dRmax, bool clusterPhotons, bool trackPhotons,
  //                  double masstarget) {
  //   vector<pair<double, double> > etaRanges;
  //   etaRanges += std::make_pair(etaMin, etaMax);
  //   FinalState inputfs;
  //   _init(inputfs, etaRanges, pTmin, pid, minmass, maxmass, dRmax, clusterPhotons, trackPhotons, masstarget);
  // }


  // ZFinder::ZFinder(const std::vector<std::pair<double, double> >& etaRanges,
  //                  double pTmin,
  //                  PdgId pid,
  //                  double minmass, const double maxmass,
  //                  double dRmax, bool clusterPhotons, bool trackPhotons,
  //                  double masstarget) {
  //   FinalState inputfs;
  //   _init(inputfs, etaRanges, pTmin, pid, minmass, maxmass, dRmax, clusterPhotons, trackPhotons, masstarget);
  // }

  void ZFinder::_init(const FinalState& inputfs, Cut fsCut,
		      // const std::vector<std::pair<double, double> >& etaRanges,
                      // double pTmin,  
		      PdgId pid,
                      double minmass, double maxmass,
                      double dRmax, bool clusterPhotons, bool trackPhotons,
                      double masstarget)
  {
    setName("ZFinder");

    _minmass = minmass;
    _maxmass = maxmass;
    _masstarget = masstarget;
    _pid = pid;
    _trackPhotons = trackPhotons;

    IdentifiedFinalState bareleptons(inputfs);
    bareleptons.acceptIdPair(pid);
    LeptonClusters leptons(inputfs, bareleptons, dRmax,
                           clusterPhotons, fsCut);
    //                           etaRanges, pTmin);
    addProjection(leptons, "LeptonClusters");

    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(*this);
    addProjection(remainingFS, "RFS");
  }


  /////////////////////////////////////////////////////


  const FinalState& ZFinder::remainingFinalState() const
  {
    return getProjection<FinalState>("RFS");
  }


  int ZFinder::compare(const Projection& p) const {
    PCmp LCcmp = mkNamedPCmp(p, "LeptonClusters");
    if (LCcmp != EQUIVALENT) return LCcmp;

    const ZFinder& other = dynamic_cast<const ZFinder&>(p);
    return (cmp(_minmass, other._minmass) || cmp(_maxmass, other._maxmass) ||
            cmp(_pid, other._pid) || cmp(_trackPhotons, other._trackPhotons));
  }


  void ZFinder::project(const Event& e) {
    clear();

    const LeptonClusters& leptons = applyProjection<LeptonClusters>(e, "LeptonClusters");

    InvMassFinalState imfs(std::make_pair(_pid, -_pid), _minmass, _maxmass, _masstarget);
    Particles tmp;
    tmp.insert(tmp.end(), leptons.clusteredLeptons().begin(), leptons.clusteredLeptons().end());
    imfs.calc(tmp);

    if (imfs.particlePairs().size() < 1) return;
    ParticlePair Zconstituents(imfs.particlePairs()[0]);
    Particle l1(Zconstituents.first), l2(Zconstituents.second);
    if (PID::threeCharge(l1)>0.0) {
      _constituents += l1, l2;
    }
    else {
      _constituents += l2, l1;
    }
    FourMomentum pZ = l1.momentum() + l2.momentum();
    const int z3charge = PID::threeCharge(l1.pdgId()) + PID::threeCharge(l2.pdgId());
    assert(z3charge == 0);

    stringstream msg;
    msg << "Z reconstructed from: \n"
        << "   " << l1.momentum() << " " << l1.pdgId() << "\n"
        << " + " << l2.momentum() << " " << l2.pdgId();
    MSG_DEBUG(msg.str());
    _bosons.push_back(Particle(PID::ZBOSON, pZ));

    // Find the LeptonClusters which survived the IMFS cut such that we can
    // extract their original particles
    foreach (const Particle& p, _constituents) {
      foreach (const ClusteredLepton& l, leptons.clusteredLeptons()) {
        if (p.pdgId()==l.pdgId() && p.momentum()==l.momentum()) {
          _theParticles.push_back(l.constituentLepton());
          if (_trackPhotons) {
            _theParticles.insert(_theParticles.end(),
                                 l.constituentPhotons().begin(), l.constituentPhotons().end());
          }
        }
      }
    }
  }


}
