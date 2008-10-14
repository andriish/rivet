// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/PDG_Hadron_Multiplicities_Ratios.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"

namespace Rivet {


  void PDG_HADRON_MULTIPLICITIES_RATIOS::analyze(const Event& e) {
    // First, veto on leptonic events by requiring at least 4 charged FS particles
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
    if (numParticles < 2) {
      getLog() << Log::DEBUG << "Failed leptonic event cut" << endl;
      vetoEvent(e);
    }
    getLog() << Log::DEBUG << "Passed leptonic event cut" << endl;

    // Get event weight for histo filling
    const double weight = e.weight();

    // Get beams and average beam momentum
    const ParticlePair& beams = applyProjection<Beam>(e, "Beams").getBeams();
    const double meanBeamMom = ( beams.first.momentum().vector3().mod() + 
                                 beams.second.momentum().vector3().mod() ) / 2.0;
    getLog() << Log::DEBUG << "Avg beam momentum = " << meanBeamMom << endl;

    // Final state of unstable particles to get particle spectra
    const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");


    if (2*meanBeamMom >= 9.5 && 2*meanBeamMom <= 10.5) {
      for (ParticleVector::const_iterator p = ufs.particles().begin(); p != ufs.particles().end(); ++p) {
        int id = abs(p->getPdgId());
        switch (id) {
           case 211:
              _weightedTotalNumPiPlus10 += weight;
              break;
           case 111:
              _hist10MeanMultiPi0->fill(_hist10MeanMultiPi0->binMean(0), weight);
              break;
           case 321:
              _hist10MeanMultiKPlus->fill(_hist10MeanMultiKPlus->binMean(0), weight);
              break;
           case 130:
           case 310:
              _hist10MeanMultiK0->fill(_hist10MeanMultiK0->binMean(0), weight);
              break;
           case 221:
              _hist10MeanMultiEta->fill(_hist10MeanMultiEta->binMean(0), weight);
              break;
           case 331:
              _hist10MeanMultiEtaPrime->fill(_hist10MeanMultiEtaPrime->binMean(0), weight);
              break;
           case 411:
              _hist10MeanMultiDPlus->fill(_hist10MeanMultiDPlus->binMean(0), weight);
              break;
           case 421:
              _hist10MeanMultiD0->fill(_hist10MeanMultiD0->binMean(0), weight);
              break;
           case 431:
              _hist10MeanMultiDPlus_s->fill(_hist10MeanMultiDPlus_s->binMean(0), weight);
              break;
           case 9010221:
              _hist10MeanMultiF0_980->fill(_hist10MeanMultiF0_980->binMean(0), weight);
              break;
           case 113:
              _hist10MeanMultiRho770_0->fill(_hist10MeanMultiRho770_0->binMean(0), weight);
              break;
           case 223:
              _hist10MeanMultiOmega782->fill(_hist10MeanMultiOmega782->binMean(0), weight);
              break;
           case 323:
              _hist10MeanMultiKStar892Plus->fill(_hist10MeanMultiKStar892Plus->binMean(0), weight);
              break;
           case 313:
              _hist10MeanMultiKStar892_0->fill(_hist10MeanMultiKStar892_0->binMean(0), weight);
              break;
           case 333:
              _hist10MeanMultiPhi1020->fill(_hist10MeanMultiPhi1020->binMean(0), weight);
              break;
           case 413:
              _hist10MeanMultiDStar2010Plus->fill(_hist10MeanMultiDStar2010Plus->binMean(0), weight);
              break;
           case 423:
              _hist10MeanMultiDStar2007_0->fill(_hist10MeanMultiDStar2007_0->binMean(0), weight);
              break;
           case 433:
              _hist10MeanMultiDStar_s2112Plus->fill(_hist10MeanMultiDStar_s2112Plus->binMean(0), weight);
              break;
           case 443:
              _hist10MeanMultiJPsi1S->fill(_hist10MeanMultiJPsi1S->binMean(0), weight);
              break;
           case 225:
              _hist10MeanMultiF2_1270->fill(_hist10MeanMultiF2_1270->binMean(0), weight);
              break;
           case 2212:
              _hist10MeanMultiP->fill(_hist10MeanMultiP->binMean(0), weight);
              break;
           case 3122:
              _hist10MeanMultiLambda->fill(_hist10MeanMultiLambda->binMean(0), weight);
              break;
           case 3212:
              _hist10MeanMultiSigma0->fill(_hist10MeanMultiSigma0->binMean(0), weight);
              break;
           case 3312:
              _hist10MeanMultiXiMinus->fill(_hist10MeanMultiXiMinus->binMean(0), weight);
              break;
           case 2224:
              _hist10MeanMultiDelta1232PlusPlus->fill(_hist10MeanMultiDelta1232PlusPlus->binMean(0), weight);
              break;
           case 3114:
              _hist10MeanMultiSigma1385Minus->fill(_hist10MeanMultiSigma1385Minus->binMean(0), weight);
              _hist10MeanMultiSigma1385PlusMinus->fill(_hist10MeanMultiSigma1385PlusMinus->binMean(0), weight);
              break;
           case 3224:
              _hist10MeanMultiSigma1385Plus->fill(_hist10MeanMultiSigma1385Plus->binMean(0), weight);
              _hist10MeanMultiSigma1385PlusMinus->fill(_hist10MeanMultiSigma1385PlusMinus->binMean(0), weight);
              break;
           case 3324:
              _hist10MeanMultiXi1530_0->fill(_hist10MeanMultiXi1530_0->binMean(0), weight);
              break;
           case 3334:
              _hist10MeanMultiOmegaMinus->fill(_hist10MeanMultiOmegaMinus->binMean(0), weight);
              break;
           case 4122:
              _hist10MeanMultiLambda_c_Plus->fill(_hist10MeanMultiLambda_c_Plus->binMean(0), weight);
              break;
           case 4222:
           case 4112:
              _hist10MeanMultiSigma_c_PlusPlus_0->fill(_hist10MeanMultiSigma_c_PlusPlus_0->binMean(0), weight);
              break;
           case 3124:
              _hist10MeanMultiLambda1520->fill(_hist10MeanMultiLambda1520->binMean(0), weight);
              break;
        }
      }
    }

    if (2*meanBeamMom >= 29 && 2*meanBeamMom <= 35) {
      for (ParticleVector::const_iterator p = ufs.particles().begin(); p != ufs.particles().end(); ++p) {
        int id = abs(p->getPdgId());
        switch (id) {
           case 211:
              _weightedTotalNumPiPlus32 += weight;
              break;
           case 111:
              _hist32MeanMultiPi0->fill(_hist32MeanMultiPi0->binMean(0), weight);
              break;
           case 321:
              _hist32MeanMultiKPlus->fill(_hist32MeanMultiKPlus->binMean(0), weight);
              break;
           case 130:
           case 310:
              _hist32MeanMultiK0->fill(_hist32MeanMultiK0->binMean(0), weight);
              break;
           case 221:
              _hist32MeanMultiEta->fill(_hist32MeanMultiEta->binMean(0), weight);
              break;
           case 331:
              _hist32MeanMultiEtaPrime->fill(_hist32MeanMultiEtaPrime->binMean(0), weight);
              break;
           case 411:
              _hist32MeanMultiDPlus->fill(_hist32MeanMultiDPlus->binMean(0), weight);
              break;
           case 421:
              _hist32MeanMultiD0->fill(_hist32MeanMultiD0->binMean(0), weight);
              break;
           case 431:
              _hist32MeanMultiDPlus_s->fill(_hist32MeanMultiDPlus_s->binMean(0), weight);
              break;
           case 9010221:
              _hist32MeanMultiF0_980->fill(_hist32MeanMultiF0_980->binMean(0), weight);
              break;
           case 113:
              _hist32MeanMultiRho770_0->fill(_hist32MeanMultiRho770_0->binMean(0), weight);
              break;
           case 323:
              _hist32MeanMultiKStar892Plus->fill(_hist32MeanMultiKStar892Plus->binMean(0), weight);
              break;
           case 313:
              _hist32MeanMultiKStar892_0->fill(_hist32MeanMultiKStar892_0->binMean(0), weight);
              break;
           case 333:
              _hist32MeanMultiPhi1020->fill(_hist32MeanMultiPhi1020->binMean(0), weight);
              break;
           case 413:
              _hist32MeanMultiDStar2010Plus->fill(_hist32MeanMultiDStar2010Plus->binMean(0), weight);
              break;
           case 423:
              _hist32MeanMultiDStar2007_0->fill(_hist32MeanMultiDStar2007_0->binMean(0), weight);
              break;
           case 225:
              _hist32MeanMultiF2_1270->fill(_hist32MeanMultiF2_1270->binMean(0), weight);
              break;
           case 325:
              _hist32MeanMultiK2Star1430Plus->fill(_hist32MeanMultiK2Star1430Plus->binMean(0), weight);
              break;
           case 315:
              _hist32MeanMultiK2Star1430_0->fill(_hist32MeanMultiK2Star1430_0->binMean(0), weight);
              break;
           case 2212:
              _hist32MeanMultiP->fill(_hist32MeanMultiP->binMean(0), weight);
              break;
           case 3122:
              _hist32MeanMultiLambda->fill(_hist32MeanMultiLambda->binMean(0), weight);
              break;
           case 3312:
              _hist32MeanMultiXiMinus->fill(_hist32MeanMultiXiMinus->binMean(0), weight);
              break;
           case 3114:
              _hist32MeanMultiSigma1385Minus->fill(_hist32MeanMultiSigma1385Minus->binMean(0), weight);
              _hist32MeanMultiSigma1385PlusMinus->fill(_hist32MeanMultiSigma1385PlusMinus->binMean(0), weight);
              break;
           case 3224:
              _hist32MeanMultiSigma1385Plus->fill(_hist32MeanMultiSigma1385Plus->binMean(0), weight);
              _hist32MeanMultiSigma1385PlusMinus->fill(_hist32MeanMultiSigma1385PlusMinus->binMean(0), weight);
              break;
           case 3334:
              _hist32MeanMultiOmegaMinus->fill(_hist32MeanMultiOmegaMinus->binMean(0), weight);
              break;
           case 4122:
              _hist32MeanMultiLambda_c_Plus->fill(_hist32MeanMultiLambda_c_Plus->binMean(0), weight);
              break;
        }
      }
    }



    if (2*meanBeamMom >= 89.5 && 2*meanBeamMom <= 91.8) {
      for (ParticleVector::const_iterator p = ufs.particles().begin(); p != ufs.particles().end(); ++p) {
        int id = abs(p->getPdgId());
        switch (id) {
           case 211:
              _weightedTotalNumPiPlus91 += weight;
              break;
           case 111:
              _hist91MeanMultiPi0->fill(_hist91MeanMultiPi0->binMean(0), weight);
              break;
           case 321:
              _hist91MeanMultiKPlus->fill(_hist91MeanMultiKPlus->binMean(0), weight);
              break;
           case 130:
           case 310:
              _hist91MeanMultiK0->fill(_hist91MeanMultiK0->binMean(0), weight);
              break;
           case 221:
              _hist91MeanMultiEta->fill(_hist91MeanMultiEta->binMean(0), weight);
              break;
           case 331:
              _hist91MeanMultiEtaPrime->fill(_hist91MeanMultiEtaPrime->binMean(0), weight);
              break;
           case 411:
              _hist91MeanMultiDPlus->fill(_hist91MeanMultiDPlus->binMean(0), weight);
              break;
           case 421:
              _hist91MeanMultiD0->fill(_hist91MeanMultiD0->binMean(0), weight);
              break;
           case 431:
              _hist91MeanMultiDPlus_s->fill(_hist91MeanMultiDPlus_s->binMean(0), weight);
              break;
           case 511:
              _hist91MeanMultiBPlus_B0_d->fill(_hist91MeanMultiBPlus_B0_d->binMean(0), weight);
              break;
           case 521:
              _hist91MeanMultiBPlus_B0_d->fill(_hist91MeanMultiBPlus_B0_d->binMean(0), weight);
              _hist91MeanMultiBPlus_u->fill(_hist91MeanMultiBPlus_u->binMean(0), weight);
              break;
           case 531:
              _hist91MeanMultiB0_s->fill(_hist91MeanMultiB0_s->binMean(0), weight);
              break;
           case 9010221:
              _hist91MeanMultiF0_980->fill(_hist91MeanMultiF0_980->binMean(0), weight);
              break;
           case 9000211:
              _hist91MeanMultiA0_980Plus->fill(_hist91MeanMultiA0_980Plus->binMean(0), weight);
              break;
           case 113:
              _hist91MeanMultiRho770_0->fill(_hist91MeanMultiRho770_0->binMean(0), weight);
              break;
           case 213:
              _hist91MeanMultiRho770Plus->fill(_hist91MeanMultiRho770Plus->binMean(0), weight);
              break;
           case 223:
              _hist91MeanMultiOmega782->fill(_hist91MeanMultiOmega782->binMean(0), weight);
              break;
           case 323:
              _hist91MeanMultiKStar892Plus->fill(_hist91MeanMultiKStar892Plus->binMean(0), weight);
              break;
           case 313:
              _hist91MeanMultiKStar892_0->fill(_hist91MeanMultiKStar892_0->binMean(0), weight);
              break;
           case 333:
              _hist91MeanMultiPhi1020->fill(_hist91MeanMultiPhi1020->binMean(0), weight);
              break;
           case 413:
              _hist91MeanMultiDStar2010Plus->fill(_hist91MeanMultiDStar2010Plus->binMean(0), weight);
              break;
           case 433:
              _hist91MeanMultiDStar_s2112Plus->fill(_hist91MeanMultiDStar_s2112Plus->binMean(0), weight);
              break;
           case 513:
           case 523:
           case 533:
              _hist91MeanMultiBStar->fill(_hist91MeanMultiBStar->binMean(0), weight);
              break;
           case 443:
              _hist91MeanMultiJPsi1S->fill(_hist91MeanMultiJPsi1S->binMean(0), weight);
              break;
           case 100443:
              _hist91MeanMultiPsi2S->fill(_hist91MeanMultiPsi2S->binMean(0), weight);
              break;
           case 553:
              _hist91MeanMultiUpsilon1S->fill(_hist91MeanMultiUpsilon1S->binMean(0), weight);
              break;
           case 20223:
              _hist91MeanMultiF1_1285->fill(_hist91MeanMultiF1_1285->binMean(0), weight);
              break;
           case 20333:
              _hist91MeanMultiF1_1420->fill(_hist91MeanMultiF1_1420->binMean(0), weight);
              break;
           case 445:
              _hist91MeanMultiChi_c1_3510->fill(_hist91MeanMultiChi_c1_3510->binMean(0), weight);
              break;
           case 225:
              _hist91MeanMultiF2_1270->fill(_hist91MeanMultiF2_1270->binMean(0), weight);
              break;
           case 335:
              _hist91MeanMultiF2Prime1525->fill(_hist91MeanMultiF2Prime1525->binMean(0), weight);
              break;
           case 315:
              _hist91MeanMultiK2Star1430_0->fill(_hist91MeanMultiK2Star1430_0->binMean(0), weight);
              break;
           case 515:
           case 525:
           case 535:
              _hist91MeanMultiBStarStar->fill(_hist91MeanMultiBStarStar->binMean(0), weight);
              break;
           case 10433:
           case 20433:
              _hist91MeanMultiDs1Plus->fill(_hist91MeanMultiDs1Plus->binMean(0), weight);
              break;
           case 435:
              _hist91MeanMultiDs2Plus->fill(_hist91MeanMultiDs2Plus->binMean(0), weight);
              break;
           case 2212:
              _hist91MeanMultiP->fill(_hist91MeanMultiP->binMean(0), weight);
              break;
           case 3122:
              _hist91MeanMultiLambda->fill(_hist91MeanMultiLambda->binMean(0), weight);
              break;
           case 3212:
              _hist91MeanMultiSigma0->fill(_hist91MeanMultiSigma0->binMean(0), weight);
              break;
           case 3112:
              _hist91MeanMultiSigmaMinus->fill(_hist91MeanMultiSigmaMinus->binMean(0), weight);
              _hist91MeanMultiSigmaPlusMinus->fill(_hist91MeanMultiSigmaPlusMinus->binMean(0), weight);
              break;
           case 3222:
              _hist91MeanMultiSigmaPlus->fill(_hist91MeanMultiSigmaPlus->binMean(0), weight);
              _hist91MeanMultiSigmaPlusMinus->fill(_hist91MeanMultiSigmaPlusMinus->binMean(0), weight);
              break;
           case 3312:
              _hist91MeanMultiXiMinus->fill(_hist91MeanMultiXiMinus->binMean(0), weight);
              break;
           case 2224:
              _hist91MeanMultiDelta1232PlusPlus->fill(_hist91MeanMultiDelta1232PlusPlus->binMean(0), weight);
              break;
           case 3114:
              _hist91MeanMultiSigma1385Minus->fill(_hist91MeanMultiSigma1385Minus->binMean(0), weight);
              _hist91MeanMultiSigma1385PlusMinus->fill(_hist91MeanMultiSigma1385PlusMinus->binMean(0), weight);
              break;
           case 3224:
              _hist91MeanMultiSigma1385Plus->fill(_hist91MeanMultiSigma1385Plus->binMean(0), weight);
              _hist91MeanMultiSigma1385PlusMinus->fill(_hist91MeanMultiSigma1385PlusMinus->binMean(0), weight);
              break;
           case 3324:
              _hist91MeanMultiXi1530_0->fill(_hist91MeanMultiXi1530_0->binMean(0), weight);
              break;
           case 3334:
              _hist91MeanMultiOmegaMinus->fill(_hist91MeanMultiOmegaMinus->binMean(0), weight);
              break;
           case 4122:
              _hist91MeanMultiLambda_c_Plus->fill(_hist91MeanMultiLambda_c_Plus->binMean(0), weight);
              break;
           case 5122:
              _hist91MeanMultiLambda_b_0->fill(_hist91MeanMultiLambda_b_0->binMean(0), weight);
              break;
           case 3124:
              _hist91MeanMultiLambda1520->fill(_hist91MeanMultiLambda1520->binMean(0), weight);
              break;
        }
      }
    }



    if (2*meanBeamMom >= 130 && 2*meanBeamMom <= 200) {
      for (ParticleVector::const_iterator p = ufs.particles().begin(); p != ufs.particles().end(); ++p) {
        int id = abs(p->getPdgId());
        switch (id) {
           case 211:
              _weightedTotalNumPiPlus165 += weight;
              break;
           case 321:
              _hist165MeanMultiKPlus->fill(_hist165MeanMultiKPlus->binMean(0), weight);
              break;
           case 130:
           case 310:
              _hist165MeanMultiK0->fill(_hist165MeanMultiK0->binMean(0), weight);
              break;
           case 2212:
              _hist165MeanMultiP->fill(_hist165MeanMultiP->binMean(0), weight);
              break;
           case 3122:
              _hist165MeanMultiLambda->fill(_hist165MeanMultiLambda->binMean(0), weight);
              break;
        }
      }
    }


  }



  void PDG_HADRON_MULTIPLICITIES_RATIOS::init() {
    _hist10MeanMultiPi0                = bookHistogram1D( 2, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Pi0 multiplicity");
    _hist10MeanMultiKPlus              = bookHistogram1D( 3, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean K+ multiplicity");
    _hist10MeanMultiK0                 = bookHistogram1D( 4, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean K0 multiplicity");
    _hist10MeanMultiEta                = bookHistogram1D( 5, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean eta multiplicity");
    _hist10MeanMultiEtaPrime           = bookHistogram1D( 6, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean eta'(958) multiplicity");
    _hist10MeanMultiDPlus              = bookHistogram1D( 7, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean D+ multiplicity");
    _hist10MeanMultiD0                 = bookHistogram1D( 8, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean D0 multiplicity");
    _hist10MeanMultiDPlus_s            = bookHistogram1D( 9, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean D+s multiplicity");
    _hist10MeanMultiF0_980             = bookHistogram1D(13, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean f0(980) multiplicity");
    _hist10MeanMultiRho770_0           = bookHistogram1D(15, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean rho(770)0 multiplicity");
    _hist10MeanMultiOmega782           = bookHistogram1D(17, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean omega(782) multiplicity");
    _hist10MeanMultiKStar892Plus       = bookHistogram1D(18, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean K*(892)+ multiplicity");
    _hist10MeanMultiKStar892_0         = bookHistogram1D(19, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean K*(892)0 multiplicity");
    _hist10MeanMultiPhi1020            = bookHistogram1D(20, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean phi(1020) multiplicity");
    _hist10MeanMultiDStar2010Plus      = bookHistogram1D(21, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean D*(2010)+ multiplicity");
    _hist10MeanMultiDStar2007_0        = bookHistogram1D(22, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean D*(2007)0 multiplicity");
    _hist10MeanMultiDStar_s2112Plus    = bookHistogram1D(23, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean D*s(2112)+ multiplicity");
    _hist10MeanMultiJPsi1S             = bookHistogram1D(25, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean J/psi(1S) multiplicity");
    _hist10MeanMultiF2_1270            = bookHistogram1D(31, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean f2(1270) multiplicity");
    _hist10MeanMultiP                  = bookHistogram1D(38, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean p multiplicity");
    _hist10MeanMultiLambda             = bookHistogram1D(39, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda multiplicity");
    _hist10MeanMultiSigma0             = bookHistogram1D(40, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma0 multiplicity");
    _hist10MeanMultiXiMinus            = bookHistogram1D(44, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Xi- multiplicity");
    _hist10MeanMultiDelta1232PlusPlus  = bookHistogram1D(45, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Delta(1232)++ multiplicity");
    _hist10MeanMultiSigma1385Minus     = bookHistogram1D(46, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma(1385)- multiplicity");
    _hist10MeanMultiSigma1385Plus      = bookHistogram1D(47, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma(1385)+ multiplicity");
    _hist10MeanMultiSigma1385PlusMinus = bookHistogram1D(48, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma(1385)+- multiplicity");
    _hist10MeanMultiXi1530_0           = bookHistogram1D(49, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Xi(1530)0 multiplicity");
    _hist10MeanMultiOmegaMinus         = bookHistogram1D(50, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Omega- multiplicity");
    _hist10MeanMultiLambda_c_Plus      = bookHistogram1D(51, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda\\_c+ multiplicity");
    _hist10MeanMultiSigma_c_PlusPlus_0 = bookHistogram1D(53, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma\\_c++, Sigma\\_c0 multiplicity");
    _hist10MeanMultiLambda1520         = bookHistogram1D(54, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda(1520) multiplicity");

    _hist32MeanMultiPi0                = bookHistogram1D( 2, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Pi0 multiplicity");
    _hist32MeanMultiKPlus              = bookHistogram1D( 3, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean K+ multiplicity");
    _hist32MeanMultiK0                 = bookHistogram1D( 4, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean K0 multiplicity");
    _hist32MeanMultiEta                = bookHistogram1D( 5, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean eta multiplicity");
    _hist32MeanMultiEtaPrime           = bookHistogram1D( 6, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean eta'(958) multiplicity");
    _hist32MeanMultiDPlus              = bookHistogram1D( 7, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean D+ multiplicity");
    _hist32MeanMultiD0                 = bookHistogram1D( 8, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean D0 multiplicity");
    _hist32MeanMultiDPlus_s            = bookHistogram1D( 9, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean D+s multiplicity");
    _hist32MeanMultiF0_980             = bookHistogram1D(13, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean f0(980) multiplicity");
    _hist32MeanMultiRho770_0           = bookHistogram1D(15, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean rho(770)0 multiplicity");
    _hist32MeanMultiKStar892Plus       = bookHistogram1D(18, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean K*(892)+ multiplicity");
    _hist32MeanMultiKStar892_0         = bookHistogram1D(19, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean K*(892)0 multiplicity");
    _hist32MeanMultiPhi1020            = bookHistogram1D(20, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean phi(1020) multiplicity");
    _hist32MeanMultiDStar2010Plus      = bookHistogram1D(21, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean D*(2010)+ multiplicity");
    _hist32MeanMultiDStar2007_0        = bookHistogram1D(22, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean D*(2007)0 multiplicity");
    _hist32MeanMultiF2_1270            = bookHistogram1D(31, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean f2(1270) multiplicity");
    _hist32MeanMultiK2Star1430Plus     = bookHistogram1D(33, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean K2*(1430)+ multiplicity");
    _hist32MeanMultiK2Star1430_0       = bookHistogram1D(34, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean K2*(1430)0 multiplicity");
    _hist32MeanMultiP                  = bookHistogram1D(38, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean p multiplicity");
    _hist32MeanMultiLambda             = bookHistogram1D(39, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda multiplicity");
    _hist32MeanMultiXiMinus            = bookHistogram1D(44, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Xi- multiplicity");
    _hist32MeanMultiSigma1385Minus     = bookHistogram1D(46, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma(1385)- multiplicity");
    _hist32MeanMultiSigma1385Plus      = bookHistogram1D(47, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma(1385)+ multiplicity");
    _hist32MeanMultiSigma1385PlusMinus = bookHistogram1D(48, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma(1385)+- multiplicity");
    _hist32MeanMultiOmegaMinus         = bookHistogram1D(50, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Omega- multiplicity");
    _hist32MeanMultiLambda_c_Plus      = bookHistogram1D(51, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda\\_c+ multiplicity");

    _hist91MeanMultiPi0                = bookHistogram1D( 2, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean Pi0 multiplicity");
    _hist91MeanMultiKPlus              = bookHistogram1D( 3, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean K+ multiplicity");
    _hist91MeanMultiK0                 = bookHistogram1D( 4, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean K0 multiplicity");
    _hist91MeanMultiEta                = bookHistogram1D( 5, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean eta multiplicity");
    _hist91MeanMultiEtaPrime           = bookHistogram1D( 6, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean eta'(958) multiplicity");
    _hist91MeanMultiDPlus              = bookHistogram1D( 7, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean D+ multiplicity");
    _hist91MeanMultiD0                 = bookHistogram1D( 8, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean D0 multiplicity");
    _hist91MeanMultiDPlus_s            = bookHistogram1D( 9, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean D+s multiplicity");
    _hist91MeanMultiBPlus_B0_d         = bookHistogram1D(10, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean B+, B0d multiplicity");
    _hist91MeanMultiBPlus_u            = bookHistogram1D(11, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean B+u multiplicity");
    _hist91MeanMultiB0_s               = bookHistogram1D(12, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean B0s multiplicity");
    _hist91MeanMultiF0_980             = bookHistogram1D(13, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean f0(980) multiplicity");
    _hist91MeanMultiA0_980Plus         = bookHistogram1D(14, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean a0(980)+ multiplicity");
    _hist91MeanMultiRho770_0           = bookHistogram1D(15, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean rho(770)0 multiplicity");
    _hist91MeanMultiRho770Plus         = bookHistogram1D(16, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean rho(770)+ multiplicity");
    _hist91MeanMultiOmega782           = bookHistogram1D(17, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean omega(782) multiplicity");
    _hist91MeanMultiKStar892Plus       = bookHistogram1D(18, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean K*(892)+ multiplicity");
    _hist91MeanMultiKStar892_0         = bookHistogram1D(19, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean K*(892)0 multiplicity");
    _hist91MeanMultiPhi1020            = bookHistogram1D(20, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean phi(1020) multiplicity");
    _hist91MeanMultiDStar2010Plus      = bookHistogram1D(21, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean D*(2010)+ multiplicity");
    _hist91MeanMultiDStar_s2112Plus    = bookHistogram1D(23, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean D*s(2112)+ multiplicity");
    _hist91MeanMultiBStar              = bookHistogram1D(24, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean B* multiplicity");
    _hist91MeanMultiJPsi1S             = bookHistogram1D(25, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean J/psi(1S) multiplicity");
    _hist91MeanMultiPsi2S              = bookHistogram1D(26, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean psi(2S) multiplicity");
    _hist91MeanMultiUpsilon1S          = bookHistogram1D(27, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Upsilon(1S) multiplicity");
    _hist91MeanMultiF1_1285            = bookHistogram1D(28, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean f1(1285) multiplicity");
    _hist91MeanMultiF1_1420            = bookHistogram1D(29, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean f1(1420) multiplicity");
    _hist91MeanMultiChi_c1_3510        = bookHistogram1D(30, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean chi\\_c1(3510) multiplicity");
    _hist91MeanMultiF2_1270            = bookHistogram1D(31, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean f2(1270) multiplicity");
    _hist91MeanMultiF2Prime1525        = bookHistogram1D(32, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean f2'(1525) multiplicity");
    _hist91MeanMultiK2Star1430_0       = bookHistogram1D(34, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean K2*(1430)0 multiplicity");
    _hist91MeanMultiBStarStar          = bookHistogram1D(35, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean B** multiplicity");
    _hist91MeanMultiDs1Plus            = bookHistogram1D(36, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Ds1+ multiplicity");
    _hist91MeanMultiDs2Plus            = bookHistogram1D(37, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Ds2+ multiplicity");
    _hist91MeanMultiP                  = bookHistogram1D(38, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean p multiplicity");
    _hist91MeanMultiLambda             = bookHistogram1D(39, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda multiplicity");
    _hist91MeanMultiSigma0             = bookHistogram1D(40, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma0 multiplicity");
    _hist91MeanMultiSigmaMinus         = bookHistogram1D(41, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma- multiplicity");
    _hist91MeanMultiSigmaPlus          = bookHistogram1D(42, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma+ multiplicity");
    _hist91MeanMultiSigmaPlusMinus     = bookHistogram1D(43, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma+- multiplicity");
    _hist91MeanMultiXiMinus            = bookHistogram1D(44, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean Xi- multiplicity");
    _hist91MeanMultiDelta1232PlusPlus  = bookHistogram1D(45, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Delta(1232)++ multiplicity");
    _hist91MeanMultiSigma1385Minus     = bookHistogram1D(46, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma(1385)- multiplicity");
    _hist91MeanMultiSigma1385Plus      = bookHistogram1D(47, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma(1385)+ multiplicity");
    _hist91MeanMultiSigma1385PlusMinus = bookHistogram1D(48, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean Sigma(1385)+- multiplicity");
    _hist91MeanMultiXi1530_0           = bookHistogram1D(49, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Xi(1530)0 multiplicity");
    _hist91MeanMultiOmegaMinus         = bookHistogram1D(50, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean Omega- multiplicity");
    _hist91MeanMultiLambda_c_Plus      = bookHistogram1D(51, 1, 3, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda\\_c+ multiplicity");
    _hist91MeanMultiLambda_b_0         = bookHistogram1D(52, 1, 1, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda\\_b0 multiplicity");
    _hist91MeanMultiLambda1520         = bookHistogram1D(54, 1, 2, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda(1520) multiplicity");

    _hist165MeanMultiKPlus             = bookHistogram1D( 3, 1, 4, "Ratio (w.r.t. Pi+/Pi-) of mean K+ multiplicity");
    _hist165MeanMultiK0                = bookHistogram1D( 4, 1, 4, "Ratio (w.r.t. Pi+/Pi-) of mean K0 multiplicity");
    _hist165MeanMultiP                 = bookHistogram1D(38, 1, 4, "Ratio (w.r.t. Pi+/Pi-) of mean p multiplicity");
    _hist165MeanMultiLambda            = bookHistogram1D(39, 1, 4, "Ratio (w.r.t. Pi+/Pi-) of mean Lambda multiplicity");
  }

  // Finalize
  void PDG_HADRON_MULTIPLICITIES_RATIOS::finalize() {
    scale(_hist10MeanMultiPi0               , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiKPlus             , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiK0                , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiEta               , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiEtaPrime          , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiDPlus             , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiD0                , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiDPlus_s           , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiF0_980            , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiRho770_0          , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiOmega782          , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiKStar892Plus      , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiKStar892_0        , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiPhi1020           , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiDStar2010Plus     , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiDStar2007_0       , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiDStar_s2112Plus   , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiJPsi1S            , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiF2_1270           , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiP                 , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiLambda            , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiSigma0            , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiXiMinus           , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiDelta1232PlusPlus , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiSigma1385Minus    , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiSigma1385Plus     , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiSigma1385PlusMinus, 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiXi1530_0          , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiOmegaMinus        , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiLambda_c_Plus     , 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiSigma_c_PlusPlus_0, 1.0/_weightedTotalNumPiPlus10);
    scale(_hist10MeanMultiLambda1520        , 1.0/_weightedTotalNumPiPlus10);

    scale(_hist32MeanMultiPi0               , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiKPlus             , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiK0                , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiEta               , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiEtaPrime          , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiDPlus             , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiD0                , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiDPlus_s           , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiF0_980            , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiRho770_0          , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiKStar892Plus      , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiKStar892_0        , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiPhi1020           , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiDStar2010Plus     , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiDStar2007_0       , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiF2_1270           , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiK2Star1430Plus    , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiK2Star1430_0      , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiP                 , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiLambda            , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiXiMinus           , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiSigma1385Minus    , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiSigma1385Plus     , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiSigma1385PlusMinus, 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiOmegaMinus        , 1.0/_weightedTotalNumPiPlus32);
    scale(_hist32MeanMultiLambda_c_Plus     , 1.0/_weightedTotalNumPiPlus32);

    scale(_hist91MeanMultiPi0               , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiKPlus             , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiK0                , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiEta               , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiEtaPrime          , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiDPlus             , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiD0                , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiDPlus_s           , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiBPlus_B0_d        , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiBPlus_u           , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiB0_s              , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiF0_980            , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiA0_980Plus        , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiRho770_0          , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiRho770Plus        , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiOmega782          , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiKStar892Plus      , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiKStar892_0        , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiPhi1020           , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiDStar2010Plus     , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiDStar_s2112Plus   , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiBStar             , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiJPsi1S            , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiPsi2S             , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiUpsilon1S         , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiF1_1285           , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiF1_1420           , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiChi_c1_3510       , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiF2_1270           , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiF2Prime1525       , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiK2Star1430_0      , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiBStarStar         , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiDs1Plus           , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiDs2Plus           , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiP                 , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiLambda            , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiSigma0            , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiSigmaMinus        , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiSigmaPlus         , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiSigmaPlusMinus    , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiXiMinus           , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiDelta1232PlusPlus , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiSigma1385Minus    , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiSigma1385Plus     , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiSigma1385PlusMinus, 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiXi1530_0          , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiOmegaMinus        , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiLambda_c_Plus     , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiLambda_b_0        , 1.0/_weightedTotalNumPiPlus91);
    scale(_hist91MeanMultiLambda1520        , 1.0/_weightedTotalNumPiPlus91);

    scale(_hist165MeanMultiKPlus            , 1.0/_weightedTotalNumPiPlus165);
    scale(_hist165MeanMultiK0               , 1.0/_weightedTotalNumPiPlus165);
    scale(_hist165MeanMultiP                , 1.0/_weightedTotalNumPiPlus165);
    scale(_hist165MeanMultiLambda           , 1.0/_weightedTotalNumPiPlus165);
  }

}
