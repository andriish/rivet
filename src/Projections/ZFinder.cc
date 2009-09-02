// -*- C++ -*-
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/ClusteredPhotons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  ZFinder::ZFinder(const FinalState& fs,
                   PdgId pid,
                   double m2_min, double m2_max,
                   double dRmax) {
    _init(fs, pid, m2_min, m2_max, dRmax);
  }


  ZFinder::ZFinder(double etaMin, double etaMax,
                   double pTmin,
                   PdgId pid,
                   double m2_min, double m2_max,
                   double dRmax) {
    vector<pair<double, double> > etaRanges;
    etaRanges += std::make_pair(etaMin, etaMax);
    _init(etaRanges, pTmin, pid, m2_min, m2_max, dRmax);
  }


  ZFinder::ZFinder(const std::vector<std::pair<double, double> >& etaRanges,
                   double pTmin,
                   PdgId pid,
                   double m2_min, const double m2_max,
                   double dRmax) {
    _init(etaRanges, pTmin, pid, m2_min, m2_max, dRmax);
  }


  void ZFinder::_init(const std::vector<std::pair<double, double> >& etaRanges,
                      double pTmin,  PdgId pid,
                      double m2_min, double m2_max,
                      double dRmax) {
    FinalState fs(etaRanges, pTmin);
    _init(fs, pid, m2_min, m2_max, dRmax);
  }


  void ZFinder::_init(const FinalState& fs,
                      PdgId pid,
                      double m2_min, double m2_max,
                      double dRmax)
  {
    setName("ZFinder");

    addProjection(fs, "FS");

    InvMassFinalState imfs(fs, std::make_pair(pid, -pid), m2_min, m2_max);
    addProjection(imfs, "IMFS");
    
    ClusteredPhotons cphotons(FinalState(), imfs, dRmax);
    addProjection(cphotons, "CPhotons");

    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(imfs);
    remainingFS.addVetoOnThisFinalState(cphotons);
    addProjection(remainingFS, "RFS");
  }


  /////////////////////////////////////////////////////


  const FinalState& ZFinder::remainingFinalState() const
  {
    return getProjection<FinalState>("RFS");
  }


  const FinalState& ZFinder::constituentsFinalState() const
  {
    return getProjection<FinalState>("IMFS");
  }

  int ZFinder::compare(const Projection& p) const {
    PCmp cmp = mkNamedPCmp(p, "IMFS");
    if (cmp != PCmp::EQUIVALENT) return cmp;

    cmp = mkNamedPCmp(p, "CPhotons");
    if (cmp != PCmp::EQUIVALENT) return cmp;

    return PCmp::EQUIVALENT;
  } 
  

  void ZFinder::project(const Event& e) {
    _theParticles.clear();

    const FinalState& imfs=applyProjection<FinalState>(e, "IMFS");
    if (imfs.particles().size() != 2) return;
    FourMomentum pZ = imfs.particles()[0].momentum() + imfs.particles()[1].momentum();
    const int z3charge = PID::threeCharge(imfs.particles()[0].pdgId()) + PID::threeCharge(imfs.particles()[1].pdgId());
    assert(z3charge == 0);

    stringstream msg;
    msg << "Z reconstructed from: " << endl
        << "   " << imfs.particles()[0].momentum() << " " << imfs.particles()[0].pdgId() << endl
        << " + " << imfs.particles()[1].momentum() << " " << imfs.particles()[1].pdgId() << endl;

    // Add in clustered photons
    const FinalState& photons = applyProjection<FinalState>(e, "CPhotons");
    foreach (const Particle& photon, photons.particles()) {
      msg << " + " << photon.momentum() << " " << photon.pdgId() << endl;
      pZ += photon.momentum();
    }
    msg << " = " << pZ;
    getLog() << Log::DEBUG << msg.str() << endl;

    Particle Z;
    Z.setMomentum(pZ);
    _theParticles.push_back(Z);
    getLog() << Log::DEBUG << name() << " found " << _theParticles.size()
             << " particles." << endl;
  }
 
 
}
