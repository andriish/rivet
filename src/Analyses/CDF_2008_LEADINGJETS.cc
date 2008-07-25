// -*- C++ -*-

// Field & Stuart underlying event analysis at CDF.
// Phys.Rev.D65:092002,2002 - no arXiv code.
// FNAL-PUB 01/211-E

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_LEADINGJETS.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {


  // Book histograms
  void CDF_2008_LEADINGJETS::init() {
    _hist_onchg   = bookProfile1D( 1, 1, 1, "Overall Number of Charged Particles");
    _hist_ocptsum = bookProfile1D( 2, 1, 1, "Overall Charged pT Sum");
    _hist_oetsum  = bookProfile1D( 3, 1, 1, "Overall ET Sum");

    //_hist_tnchg   = bookProfile1D( 1, 1, 1, "Toward Region Charged Particle Density");
    //_hist_pnchg   = bookProfile1D( 2, 1, 1, "Transverse Region Charged Particle Density");
    //_hist_anchg   = bookProfile1D( 3, 1, 1, "Away Region Charged Particle Density");
    //_hist_tcptsum = bookProfile1D( 4, 1, 1, "Toward Region Charged pT Sum Density");
    //_hist_pcptsum = bookProfile1D( 5, 1, 1, "Transverse Region Charged pT Sum Density");
    //_hist_acptsum = bookProfile1D( 6, 1, 1, "Away Region Charged pT Sum Density");
    //_hist_tcptave = bookProfile1D( 7, 1, 1, "Toward Region Charged pT Average");
    //_hist_pcptave = bookProfile1D( 8, 1, 1, "Transverse Region Charged pT Average");
    //_hist_acptave = bookProfile1D( 9, 1, 1, "Away Region Charged pT Average");
    //_hist_tcptmax = bookProfile1D(10, 1, 1, "Toward Region Charged pT Maximum");
    //_hist_pcptmax = bookProfile1D(11, 1, 1, "Transverse Region Charged pT Maximum");
    //_hist_acptmax = bookProfile1D(12, 1, 1, "Away Region Charged pT Maximum");
  }


  // Do the analysis
  void CDF_2008_LEADINGJETS::analyze(const Event& e) {
    Log log = getLog();

    const FinalState& fsj = applyProjection<FinalState>(e, "FSJ");
    if (fsj.particles().size() < 1) {
      getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
      vetoEvent(e);
    }

    const FastJets& jetpro = applyProjection<FastJets>(e, "MidpointJets");
    const PseudoJets& jets = jetpro.getPseudoJetsByPt();

    getLog() << Log::DEBUG << "jet multiplicity = " << jets.size() << endl;

    // We require the leading jet to be within |eta|<2
    if (jets.size() < 1 || fabs(jets[0].eta()) >= 2) {
      getLog() << Log::DEBUG << "Failed jet cut" << endl;
      vetoEvent(e);
    }

    const double jetphi = jets[0].phi();
    const double jetpT  = jets[0].perp();

    getLog() << Log::DEBUG << "leading jet: pT = " << jetpT
      << ", eta = " << jets[0].eta()
      << ", phi = " << jetphi << endl;

    // Get the event weight
    const double weight = e.weight();

    // Get the final states to work with for filling the distributions
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(e, "CFS");

//    getLog() << Log::DEBUG << fsj.particles().size() << "   "
//                           <<  fs.particles().size() << "   "
//                           << cfs.particles().size() << endl;

    size_t   numOverall(0),     numToward(0),     numTrans1(0),     numTrans2(0),     numAway(0)  ;
    double ptSumOverall(0.0), ptSumToward(0.0), ptSumTrans1(0.0), ptSumTrans2(0.0), ptSumAway(0.0);
    double EtSumOverall(0.0), EtSumToward(0.0), EtSumTrans1(0.0), EtSumTrans2(0.0), EtSumAway(0.0);
    double ptMaxOverall(0.0), ptMaxToward(0.0), ptMaxTrans1(0.0), ptMaxTrans2(0.0), ptMaxAway(0.0);

    // Calculate all the charged stuff
    for (ParticleVector::const_iterator p = cfs.particles().begin(); p != cfs.particles().end(); ++p) {
      const double deltaPhi = delta_phi(p->getMomentum().azimuthalAngle(), jetphi);
      const double pT = p->getMomentum().pT();
      const double phi = p->getMomentum().azimuthalAngle();

      // Jets come with phi in [0 .. 2*Pi]. Particles come in [-Pi .. Pi].
      // Lovely, isn't it?
      double rotatedphi = phi - jetphi;
      while (rotatedphi < 0) rotatedphi += 2*PI;

      ptSumOverall += pT;
      ++numOverall;
      if (pT > ptMaxOverall)
        ptMaxOverall = pT;

      if (deltaPhi < PI/3.0) {
        ptSumToward += pT;
        ++numToward;
        if (pT > ptMaxToward)
          ptMaxToward = pT;
      }
      else if (deltaPhi < 2*PI/3.0) {
        if (rotatedphi <= PI) {
          ptSumTrans1 += pT;
          ++numTrans1;
          if (pT > ptMaxTrans1)
            ptMaxTrans1 = pT;
        }
        else {
          ptSumTrans2 += pT;
          ++numTrans2;
          if (pT > ptMaxTrans2)
            ptMaxTrans2 = pT;
        }
      }
      else {
        ptSumAway += pT;
        ++numAway;
        if (pT > ptMaxAway)
          ptMaxAway = pT;
      }
    } // end charged particle loop

    // And now the same business for all particles (including neutrals)
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      const double deltaPhi = delta_phi(p->getMomentum().azimuthalAngle(), jetphi);
      const double ET = p->getMomentum().Et();
      const double phi = p->getMomentum().azimuthalAngle();

      // Jets come with phi in [0 .. 2*Pi]. Particles come in [-Pi .. Pi].
      // Lovely, isn't it?
      double rotatedphi = phi - jetphi;
      while (rotatedphi < 0) rotatedphi += 2*PI;

      EtSumOverall += ET;

      if (deltaPhi < PI/3.0) {
        EtSumToward += ET;
      }
      else if (deltaPhi < 2*PI/3.0) {
        if (rotatedphi <= PI) {
          EtSumTrans1 += ET;
        }
        else {
          EtSumTrans2 += ET;
        }
      }
      else {
        EtSumAway += ET;
      }
    } // end all particle loop

    _hist_onchg->fill(jetpT, numOverall, weight);
    _hist_ocptsum->fill(jetpT, ptSumOverall, weight);
    _hist_oetsum->fill(jetpT, EtSumOverall, weight);
////    _hist_tnchg->fill(pTZ, numToward/(4*PI/3), weight);
////    _hist_pnchg->fill(pTZ, numTrans/(4*PI/3), weight);
////    _hist_anchg->fill(pTZ, numAway/(4*PI/3), weight);
////    _hist_tcptsum->fill(pTZ, ptSumToward/(4*PI/3), weight);
////    _hist_pcptsum->fill(pTZ, ptSumTrans/(4*PI/3), weight);
////    _hist_acptsum->fill(pTZ, ptSumAway/(4*PI/3), weight);
////    if (numToward > 0) {
////      _hist_tcptave->fill(pTZ, ptSumToward/numToward, weight);
////      _hist_tcptmax->fill(pTZ, ptMaxToward, weight);
////    }
////    if (numTrans > 0) {
////      _hist_pcptave->fill(pTZ, ptSumTrans/numTrans, weight);
////      _hist_pcptmax->fill(pTZ, ptMaxTrans, weight);
////    }
////    if (numAway > 0) {
////      _hist_acptave->fill(pTZ, ptSumAway/numAway, weight);
////      _hist_acptmax->fill(pTZ, ptMaxAway, weight);
////    }
  }


  void CDF_2008_LEADINGJETS::finalize() { 
  }


}
