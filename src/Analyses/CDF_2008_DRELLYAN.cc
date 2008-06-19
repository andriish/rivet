// -*- C++ -*-

// Field & Stuart underlying event analysis at CDF.
// Phys.Rev.D65:092002,2002 - no arXiv code.
// FNAL-PUB 01/211-E

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_DRELLYAN.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {


  // Book histograms
  void CDF_2008_DRELLYAN::init() {
    _hist_tnchg   = bookProfile1D( 1, 1, 1, "Toward Region Charged Particle Density");
    _hist_pnchg   = bookProfile1D( 2, 1, 1, "Transverse Region Charged Particle Density");
    _hist_anchg   = bookProfile1D( 3, 1, 1, "Away Region Charged Particle Density");
    _hist_tcptsum = bookProfile1D( 4, 1, 1, "Toward Region Charged pT Sum Density");
    _hist_pcptsum = bookProfile1D( 5, 1, 1, "Transverse Region Charged pT Sum Density");
    _hist_acptsum = bookProfile1D( 6, 1, 1, "Away Region Charged pT Sum Density");
    _hist_tcptave = bookProfile1D( 7, 1, 1, "Toward Region Charged pT Average");
    _hist_pcptave = bookProfile1D( 8, 1, 1, "Transverse Region Charged pT Average");
    _hist_acptave = bookProfile1D( 9, 1, 1, "Away Region Charged pT Average");
    _hist_tcptmax = bookProfile1D(10, 1, 1, "Toward Region Charged pT Maximum");
    _hist_pcptmax = bookProfile1D(11, 1, 1, "Transverse Region Charged pT Maximum");
    _hist_acptmax = bookProfile1D(12, 1, 1, "Away Region Charged pT Maximum");
  }


  // Do the analysis
  void CDF_2008_DRELLYAN::analyze(const Event& e) {
    Log log = getLog();

    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
    if (numParticles < 1) {
      getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
      vetoEvent(e);
    }

    const ParticleVector& leptons = applyProjection<ChargedLeptons>(e, "CL").chargedLeptons();
    // Get the event weight
    const double weight = e.weight();

    getLog() << Log::DEBUG << "lepton multiplicity = " << leptons.size() << endl;
    if (leptons.size() != 2 || leptons[0].getPdgId() != -leptons[1].getPdgId() )
      vetoEvent(e);

    FourVector dilepton = leptons[0].getMomentum() + leptons[1].getMomentum();

    if (mass(dilepton) < 70 || mass(dilepton) > 110)
      vetoEvent(e);
    getLog() << Log::DEBUG << "dilepton mass = " << mass(dilepton) << endl;
    getLog() << Log::DEBUG << "dilepton pT   = " << pT(dilepton) << endl;

    double ptSumToward(0.0), ptSumAway(0.0), ptSumTrans(0.0);
    double ptMaxToward(0.0), ptMaxAway(0.0), ptMaxTrans(0.0);
    size_t numToward(0), numTrans(0), numAway(0);
    const double phiZ = azimuthalAngle(dilepton);
    const double pTZ  = pT(dilepton);
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      // don't use the leptons:
      if (abs(p->getPdgId()) < 20)
        continue;
      const double deltaPhi = delta_phi(p->getMomentum().azimuthalAngle(), phiZ);
      const double pT = p->getMomentum().pT();
      if (deltaPhi < PI/3.0) {
        ptSumToward += pT;
        ++numToward;
        if (pT > ptMaxToward)
          ptMaxToward = pT;
      } else if (deltaPhi < 2*PI/3.0) {
        ptSumTrans += pT;
        ++numTrans;
        if (pT > ptMaxTrans)
          ptMaxTrans = pT;
      } else {
        ptSumAway += pT;
        ++numAway;
        if (pT > ptMaxAway)
          ptMaxAway = pT;
      }

    }
    _hist_tnchg->fill(pTZ, numToward/(4*PI/3), weight);
    _hist_pnchg->fill(pTZ, numTrans/(4*PI/3), weight);
    _hist_anchg->fill(pTZ, numAway/(4*PI/3), weight);
    _hist_tcptsum->fill(pTZ, ptSumToward/(4*PI/3), weight);
    _hist_pcptsum->fill(pTZ, ptSumTrans/(4*PI/3), weight);
    _hist_acptsum->fill(pTZ, ptSumAway/(4*PI/3), weight);
    if (numToward > 0) {
      _hist_tcptave->fill(pTZ, ptSumToward/numToward, weight);
      _hist_tcptmax->fill(pTZ, ptMaxToward, weight);
    }
    if (numTrans > 0) {
      _hist_pcptave->fill(pTZ, ptSumTrans/numTrans, weight);
      _hist_pcptmax->fill(pTZ, ptMaxTrans, weight);
    }
    if (numAway > 0) {
      _hist_acptave->fill(pTZ, ptSumAway/numAway, weight);
      _hist_acptmax->fill(pTZ, ptMaxAway, weight);
    }
  }


  void CDF_2008_DRELLYAN::finalize() { 
  }


}
