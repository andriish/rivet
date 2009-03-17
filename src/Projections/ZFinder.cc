// -*- C++ -*-
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/ClusteredPhotons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  ZFinder::ZFinder(const std::vector<std::pair<double, double> >& etaRanges,
                   const double& pTmin,  const PdgId& pid,
                   const double& m2_min, const double& m2_max,
                   double dRmax)
  {
    setName("ZFinder");

    FinalState fs(etaRanges, pTmin);
    InvMassFinalState imfs(fs, make_pair(pid, -pid), m2_min, m2_max);
    addProjection(imfs, "IMFS");
    
    ClusteredPhotons cphotons(FinalState(), imfs, dRmax);
    addProjection(cphotons, "CPhotons");

    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(imfs);
    remainingFS.addVetoOnThisFinalState(cphotons);
    addProjection(remainingFS, "RFS");
  }

  const FinalState& ZFinder::remainingFinalState() const
  {
    return getProjection<FinalState>("RFS");
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
    if (imfs.particles().size()!=2) return;

    const FinalState& photons=applyProjection<FinalState>(e, "CPhotons");

    Particle Z;
    FourMomentum pZ=imfs.particles()[0].momentum()+imfs.particles()[1].momentum();
    foreach (const Particle& photon, photons.particles()) {
      pZ+=photon.momentum();
    }
    Z.setMomentum(pZ);

    _theParticles.push_back(Z);
    getLog() << Log::DEBUG << name() <<" found " << _theParticles.size()
             <<" particles." << endl;
  }
 
 
}
