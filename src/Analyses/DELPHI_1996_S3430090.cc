// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/DELPHI_1996_S3430090.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"

namespace Rivet {



  DELPHI_1996_S3430090::DELPHI_1996_S3430090() 
  {
    setBeams(ELECTRON, POSITRON); 
    addProjection(Beam(), "Beams");
    const ChargedFinalState cfs;
    addProjection(cfs, "FS");
    addProjection(UnstableFinalState(), "UFS");
    #ifdef HAVE_JADE
    addProjection(FastJets(cfs, FastJets::JADE, 0.7), "JadeJets");
    addProjection(FastJets(cfs, FastJets::DURHAM, 0.7), "DurhamJets");
    #endif
    addProjection(Sphericity(cfs), "Sphericity");
    addProjection(ParisiTensor(cfs), "Parisi");
    const Thrust thrust(cfs);
    addProjection(thrust, "Thrust");
    addProjection(Hemispheres(thrust), "Hemispheres");
    _weightedTotalPartNum = 0;
  }



  void DELPHI_1996_S3430090::analyze(const Event& e) {
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
    _weightedTotalPartNum += numParticles * weight;

    // Get beams and average beam momentum
    const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
    const double meanBeamMom = ( beams.first.momentum().vector3().mod() + 
                                 beams.second.momentum().vector3().mod() ) / 2.0;
    getLog() << Log::DEBUG << "Avg beam momentum = " << meanBeamMom << endl;

    // Thrusts
    getLog() << Log::DEBUG << "Calculating thrust" << endl;
    const Thrust& thrust = applyProjection<Thrust>(e, "Thrust");
    _hist1MinusT->fill(1 - thrust.thrust(), weight); 
    _histTMajor->fill(thrust.thrustMajor(), weight); 
    _histTMinor->fill(thrust.thrustMinor(), weight); 
    _histOblateness->fill(thrust.oblateness(), weight);

    // Jets
    #ifdef HAVE_JADE
    getLog() << Log::DEBUG << "Using FastJet JADE patch to make diff jet rate plots:" << endl;
    const FastJets& durjet = applyProjection<FastJets>(e, "DurhamJets");
    if (durjet.clusterSeq()) {
      _histDiffRate2Durham->fill(durjet.clusterSeq()->exclusive_dmerge(2), weight); 
      _histDiffRate3Durham->fill(durjet.clusterSeq()->exclusive_dmerge(3), weight); 
      _histDiffRate4Durham->fill(durjet.clusterSeq()->exclusive_dmerge(4), weight); 
    }
    const FastJets& jadejet = applyProjection<FastJets>(e, "JadeJets");
    if (jadejet.clusterSeq()) {
      _histDiffRate2Jade->fill(jadejet.clusterSeq()->exclusive_dmerge(2), weight); 
      _histDiffRate3Jade->fill(jadejet.clusterSeq()->exclusive_dmerge(3), weight); 
      _histDiffRate4Jade->fill(jadejet.clusterSeq()->exclusive_dmerge(4), weight); 
    }
    #endif

    // Sphericities
    getLog() << Log::DEBUG << "Calculating sphericity" << endl;
    const Sphericity& sphericity = applyProjection<Sphericity>(e, "Sphericity");
    _histSphericity->fill(sphericity.sphericity(), weight); 
    _histAplanarity->fill(sphericity.aplanarity(), weight); 
    _histPlanarity->fill(sphericity.planarity(), weight); 

    // C & D params
    getLog() << Log::DEBUG << "Calculating Parisi params" << endl;
    const ParisiTensor& parisi = applyProjection<ParisiTensor>(e, "Parisi");
    _histCParam->fill(parisi.C(), weight);
    _histDParam->fill(parisi.D(), weight);

    // Hemispheres
    getLog() << Log::DEBUG << "Calculating hemisphere variables" << endl;
    const Hemispheres& hemi = applyProjection<Hemispheres>(e, "Hemispheres");
    _histHemiMassH->fill(hemi.getScaledM2high(), weight); 
    _histHemiMassL->fill(hemi.getScaledM2low(), weight); 
    _histHemiMassD->fill(hemi.getScaledM2diff(), weight); 
    _histHemiBroadW->fill(hemi.getBmax(), weight); 
    _histHemiBroadN->fill(hemi.getBmin(), weight); 
    _histHemiBroadT->fill(hemi.getBsum(), weight); 
    _histHemiBroadD->fill(hemi.getBdiff(), weight); 

    // Iterate over all the charged final state particles.
    double Evis = 0.0;
    double Evis2 = 0.0;
    getLog() << Log::DEBUG << "About to iterate over charged FS particles" << endl;
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      // Get momentum and energy of each particle.
      const Vector3 mom3 = p->momentum().vector3();
      const double energy = p->momentum().E();
      Evis += energy;

      // Scaled momenta.
      const double mom = mom3.mod();
      const double scaledMom = mom/meanBeamMom;
      const double logInvScaledMom = -std::log(scaledMom);
      _histLogScaledMom->fill(logInvScaledMom, weight); 
      _histScaledMom->fill(scaledMom, weight); 

      // Get momenta components w.r.t. thrust and sphericity.
      const double momT = dot(thrust.thrustAxis(), mom3);
      const double momS = dot(sphericity.sphericityAxis(), mom3);
      const double pTinT = dot(mom3, thrust.thrustMajorAxis());
      const double pToutT = dot(mom3, thrust.thrustMinorAxis());
      const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
      const double pToutS = dot(mom3, sphericity.sphericityMinorAxis());
      const double pT = sqrt(pow(pTinT, 2) + pow(pToutT, 2));
      _histPtTIn->fill(fabs(pTinT/GeV), weight);
      _histPtTOut->fill(fabs(pToutT/GeV), weight);
      _histPtSIn->fill(fabs(pTinS/GeV), weight);
      _histPtSOut->fill(fabs(pToutS/GeV), weight);
      _histPtVsXp->fill(scaledMom, fabs(pT/GeV), weight);
      _histPtTOutVsXp->fill(scaledMom, fabs(pToutT/GeV), weight);

      // Calculate rapidities w.r.t. thrust and sphericity.
      const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
      const double rapidityS = 0.5 * std::log((energy + momS) / (energy - momS));
      _histRapidityT->fill(rapidityT, weight); 
      _histRapidityS->fill(rapidityS, weight); 
    }
    Evis2 = Evis*Evis;

    for (ParticleVector::const_iterator p_i = fs.particles().begin(); p_i != fs.particles().end(); ++p_i) {
      for (ParticleVector::const_iterator p_j = p_i; p_j != fs.particles().end(); ++p_j) {
        if (p_i == p_j) continue;
        const Vector3 mom3_i = p_i->momentum().vector3();
        const Vector3 mom3_j = p_j->momentum().vector3();
        const double energy_i = p_i->momentum().E();
        const double energy_j = p_j->momentum().E();
        const double cosij = dot(mom3_i.unit(), mom3_j.unit());
        const double eec = (energy_i*energy_j) / Evis2;
        _histEEC->fill(cosij, eec*weight);
        _histAEEC->fill( cosij,  eec*weight);
        _histAEEC->fill(-cosij, -eec*weight);
      }
    }

    _histMultiCharged->fill(_histMultiCharged->binMean(0), numParticles*weight);


    // Final state of unstable particles to get particle spectra
    const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");

    for (ParticleVector::const_iterator p = ufs.particles().begin(); p != ufs.particles().end(); ++p) {
      int id = abs(p->pdgId());
      switch (id) {
         case 211:
            _histMultiPiPlus->fill(_histMultiPiPlus->binMean(0), weight);
            break;
         case 111:
            _histMultiPi0->fill(_histMultiPi0->binMean(0), weight);
            break;
         case 321:
            _histMultiKPlus->fill(_histMultiKPlus->binMean(0), weight);
            break;
         case 130:
         case 310:
            _histMultiK0->fill(_histMultiK0->binMean(0), weight);
            break;
         case 221:
            _histMultiEta->fill(_histMultiEta->binMean(0), weight);
            break;
         case 331:
            _histMultiEtaPrime->fill(_histMultiEtaPrime->binMean(0), weight);
            break;
         case 411:
            _histMultiDPlus->fill(_histMultiDPlus->binMean(0), weight);
            break;
         case 421:
            _histMultiD0->fill(_histMultiD0->binMean(0), weight);
            break;
         case 511:
         case 521:
         case 531:
            _histMultiBPlus0->fill(_histMultiBPlus0->binMean(0), weight);
            break;
         case 9010221:
            _histMultiF0->fill(_histMultiF0->binMean(0), weight);
            break;
         case 113:
            _histMultiRho->fill(_histMultiRho->binMean(0), weight);
            break;
         case 323:
            _histMultiKStar892Plus->fill(_histMultiKStar892Plus->binMean(0), weight);
            break;
         case 313:
            _histMultiKStar892_0->fill(_histMultiKStar892_0->binMean(0), weight);
            break;
         case 333:
            _histMultiPhi->fill(_histMultiPhi->binMean(0), weight);
            break;
         case 413:
            _histMultiDStar2010Plus->fill(_histMultiDStar2010Plus->binMean(0), weight);
            break;
         case 225:
            _histMultiF2->fill(_histMultiF2->binMean(0), weight);
            break;
         case 315:
            _histMultiK2Star1430_0->fill(_histMultiK2Star1430_0->binMean(0), weight);
            break;
         case 2212:
            _histMultiP->fill(_histMultiP->binMean(0), weight);
            break;
         case 3122:
            _histMultiLambda0->fill(_histMultiLambda0->binMean(0), weight);
            break;
         case 3312:
            _histMultiXiMinus->fill(_histMultiXiMinus->binMean(0), weight);
            break;
         case 3334:
            _histMultiOmegaMinus->fill(_histMultiOmegaMinus->binMean(0), weight);
            break;
         case 2224:
            _histMultiDeltaPlusPlus->fill(_histMultiDeltaPlusPlus->binMean(0), weight);
            break;
         case 3114:
            _histMultiSigma1385Plus->fill(_histMultiSigma1385Plus->binMean(0), weight);
            break;
         case 3324:
            _histMultiXi1530_0->fill(_histMultiXi1530_0->binMean(0), weight);
            break;
         case 5122:
            _histMultiLambdaB0->fill(_histMultiLambdaB0->binMean(0), weight);
            break;
      }
    }

  }



  void DELPHI_1996_S3430090::init() {
    _histPtTIn       = bookHistogram1D(1, 1, 1, "In-plane $p_T$ in GeV w.r.t. thrust axes (charged)");
    _histPtTOut      = bookHistogram1D(2, 1, 1, "Out-of-plane $p_T$ in GeV w.r.t. thrust axes (charged)");
    _histPtSIn       = bookHistogram1D(3, 1, 1, "In-plane $p_T$ in GeV w.r.t. sphericity axes (charged)");
    _histPtSOut      = bookHistogram1D(4, 1, 1, "Out-of-plane $p_T$ in GeV w.r.t. sphericity axes (charged)");

    _histRapidityT   = bookHistogram1D(5, 1, 1, "Rapidity w.r.t. thrust axes, $y_T$ (charged)");
    _histRapidityS   = bookHistogram1D(6, 1, 1, "Rapidity w.r.t. sphericity axes, $y_S$ (charged)");

    _histScaledMom    = bookHistogram1D(7, 1, 1, "Scaled momentum, $x_p = |p|/|p_\\text{beam}|$ (charged)");
    _histLogScaledMom = bookHistogram1D(8, 1, 1, "Log of scaled momentum, $\\log(1/x_p)$ (charged)");

    _histPtTOutVsXp   = bookProfile1D(9,  1, 1, "Mean out-of-plane $p_T$ in GeV w.r.t. thrust axes vs. $x_p$ (charged)"); // binned in Xp
    _histPtVsXp       = bookProfile1D(10, 1, 1, "Mean $p_T$ in GeV vs. $x_p$ (charged)"); // binned in Xp

    _hist1MinusT     = bookHistogram1D(11, 1, 1, "$1-\\text{Thrust}$ (charged)");
    _histTMajor      = bookHistogram1D(12, 1, 1, "Thrust major, $M$ (charged)");
    _histTMinor      = bookHistogram1D(13, 1, 1, "Thrust minor, $m$ (charged)");
    _histOblateness  = bookHistogram1D(14, 1, 1, "Oblateness = $M - m$ (charged)");

    _histSphericity  = bookHistogram1D(15, 1, 1, "Sphericity, $S$ (charged)");
    _histAplanarity  = bookHistogram1D(16, 1, 1, "Aplanarity, $A$ (charged)");
    _histPlanarity   = bookHistogram1D(17, 1, 1, "Planarity, $P$ (charged)");

    _histCParam      = bookHistogram1D(18, 1, 1, "$C$ parameter (charged)");
    _histDParam      = bookHistogram1D(19, 1, 1, "$D$ parameter (charged)");

    _histHemiMassH   = bookHistogram1D(20, 1, 1, "Heavy hemisphere masses, $M_h^2/E_\\text{vis}^2$ (charged)");
    _histHemiMassL   = bookHistogram1D(21, 1, 1, "Light hemisphere masses, $M_l^2/E_\\text{vis}^2$ (charged)");
    _histHemiMassD   = bookHistogram1D(22, 1, 1, "Difference in hemisphere masses, $M_d^2/E_\\text{vis}^2$ (charged)");

    _histHemiBroadW  = bookHistogram1D(23, 1, 1, "Wide hemisphere broadening, $B_\\text{max}$ (charged)");
    _histHemiBroadN  = bookHistogram1D(24, 1, 1, "Narrow hemisphere broadening, $B_\\text{min}$ (charged)");
    _histHemiBroadT  = bookHistogram1D(25, 1, 1, "Total hemisphere broadening, $B_\\text{sum}$ (charged)");
    _histHemiBroadD  = bookHistogram1D(26, 1, 1, "Difference in hemisphere broadening, $B_\\text{diff}$ (charged)");

    #ifndef HAVE_JADE
    getLog() << Log::WARN << "No FastJet JADE patch, so not making any diff jet rate plots." << endl;
    #endif
    _histDiffRate2Durham  = bookHistogram1D(27, 1, 1, "Differential 2-jet rate with Durham algorithm, $D_2^\\text{Durham}$ (charged)"); // binned in y_cut
    _histDiffRate2Jade    = bookHistogram1D(28, 1, 1, "Differential 2-jet rate with Jade algorithm, $D_2^\\text{Jade}$ (charged)"); // binned in y_cut
    _histDiffRate3Durham  = bookHistogram1D(29, 1, 1, "Differential 3-jet rate with Durham algorithm, $D_3^\\text{Durham}$ (charged)"); // binned in y_cut
    _histDiffRate3Jade    = bookHistogram1D(30, 1, 1, "Differential 3-jet rate with Jade algorithm, $D_3^\\text{Jade}$ (charged)"); // binned in y_cut
    _histDiffRate4Durham  = bookHistogram1D(31, 1, 1, "Differential 4-jet rate with Durham algorithm, $D_4^\\text{Durham}$ (charged)"); // binned in y_cut
    _histDiffRate4Jade    = bookHistogram1D(32, 1, 1, "Differential 4-jet rate with Jade algorithm, $D_4^\\text{Jade}$ (charged)"); // binned in y_cut
    _histEEC              = bookHistogram1D(33, 1, 1, "Energy-energy correlation, EEC (charged)"); // binned in cos(chi)
    _histAEEC             = bookHistogram1D(34, 1, 1, "Asymmetry of the energy-energy correlation, AEEC (charged)"); // binned in cos(chi)
    _histMultiCharged     = bookHistogram1D(35, 1, 1, "Mean charged multiplicity");

    _histMultiPiPlus        = bookHistogram1D(36, 1, 1, "Mean pi+/pi- multiplicity");
    _histMultiPi0           = bookHistogram1D(36, 1, 2, "Mean pi0 multiplicity");
    _histMultiKPlus         = bookHistogram1D(36, 1, 3, "Mean K+/K- multiplicity");
    _histMultiK0            = bookHistogram1D(36, 1, 4, "Mean K0 multiplicity");
    _histMultiEta           = bookHistogram1D(36, 1, 5, "Mean eta multiplicity");
    _histMultiEtaPrime      = bookHistogram1D(36, 1, 6, "Mean etaprime multiplicity");
    _histMultiDPlus         = bookHistogram1D(36, 1, 7, "Mean D+ multiplicity");
    _histMultiD0            = bookHistogram1D(36, 1, 8, "Mean D0 multiplicity");
    _histMultiBPlus0        = bookHistogram1D(36, 1, 9, "Mean B+/B-/B0 multiplicity");

    _histMultiF0            = bookHistogram1D(37, 1, 1, "Mean f0(980) multiplicity");

    _histMultiRho           = bookHistogram1D(38, 1, 1, "Mean rho multiplicity");
    _histMultiKStar892Plus  = bookHistogram1D(38, 1, 2, "Mean K*(892)+/K*(892)- multiplicity");
    _histMultiKStar892_0    = bookHistogram1D(38, 1, 3, "Mean K*(892)0 multiplicity");
    _histMultiPhi           = bookHistogram1D(38, 1, 4, "Mean phi multiplicity");
    _histMultiDStar2010Plus = bookHistogram1D(38, 1, 5, "Mean D*(2010)+/D*(2010)- multiplicity");

    _histMultiF2            = bookHistogram1D(39, 1, 1, "Mean f2(1270) multiplicity");
    _histMultiK2Star1430_0  = bookHistogram1D(39, 1, 2, "Mean K2*(1430)0 multiplicity");

    _histMultiP             = bookHistogram1D(40, 1, 1, "Mean p multiplicity");
    _histMultiLambda0       = bookHistogram1D(40, 1, 2, "Mean Lambda0 multiplicity");
    _histMultiXiMinus       = bookHistogram1D(40, 1, 3, "Mean Xi- multiplicity");
    _histMultiOmegaMinus    = bookHistogram1D(40, 1, 4, "Mean Omega- multiplicity");
    _histMultiDeltaPlusPlus = bookHistogram1D(40, 1, 5, "Mean Delta(1232)++ multiplicity");
    _histMultiSigma1385Plus = bookHistogram1D(40, 1, 6, "Mean Sigma(1385)+/Sigma(1385)- multiplicity");
    _histMultiXi1530_0      = bookHistogram1D(40, 1, 7, "Mean Xi(1530)0 multiplicity");
    _histMultiLambdaB0      = bookHistogram1D(40, 1, 8, "Mean Lambda\\_b0 multiplicity");
  }



  // Finalize
  void DELPHI_1996_S3430090::finalize() { 
    // Normalize inclusive single particle distributions to the average number 
    // of charged particles per event.
    const double avgNumParts = _weightedTotalPartNum / sumOfWeights();

    normalize(_histPtTIn, avgNumParts);
    normalize(_histPtTOut, avgNumParts); 
    normalize(_histPtSIn, avgNumParts);
    normalize(_histPtSOut, avgNumParts); 

    normalize(_histRapidityT, avgNumParts); 
    normalize(_histRapidityS, avgNumParts); 

    normalize(_histLogScaledMom, avgNumParts);
    normalize(_histScaledMom, avgNumParts); 

    scale(_histEEC, 1.0/sumOfWeights());
    scale(_histAEEC, 1.0/sumOfWeights());
    scale(_histMultiCharged, 1.0/sumOfWeights());

    scale(_histMultiPiPlus       , 1.0/sumOfWeights());
    scale(_histMultiPi0          , 1.0/sumOfWeights());
    scale(_histMultiKPlus        , 1.0/sumOfWeights());
    scale(_histMultiK0           , 1.0/sumOfWeights());
    scale(_histMultiEta          , 1.0/sumOfWeights());
    scale(_histMultiEtaPrime     , 1.0/sumOfWeights());
    scale(_histMultiDPlus        , 1.0/sumOfWeights());
    scale(_histMultiD0           , 1.0/sumOfWeights());
    scale(_histMultiBPlus0       , 1.0/sumOfWeights());

    scale(_histMultiF0           , 1.0/sumOfWeights());

    scale(_histMultiRho          , 1.0/sumOfWeights());
    scale(_histMultiKStar892Plus , 1.0/sumOfWeights());
    scale(_histMultiKStar892_0   , 1.0/sumOfWeights());
    scale(_histMultiPhi          , 1.0/sumOfWeights());
    scale(_histMultiDStar2010Plus, 1.0/sumOfWeights());

    scale(_histMultiF2           , 1.0/sumOfWeights());
    scale(_histMultiK2Star1430_0 , 1.0/sumOfWeights());

    scale(_histMultiP             , 1.0/sumOfWeights());
    scale(_histMultiLambda0       , 1.0/sumOfWeights());
    scale(_histMultiXiMinus       , 1.0/sumOfWeights());
    scale(_histMultiOmegaMinus    , 1.0/sumOfWeights());
    scale(_histMultiDeltaPlusPlus , 1.0/sumOfWeights());
    scale(_histMultiSigma1385Plus , 1.0/sumOfWeights());
    scale(_histMultiXi1530_0      , 1.0/sumOfWeights());
    scale(_histMultiLambdaB0      , 1.0/sumOfWeights());

    normalize(_hist1MinusT); 
    normalize(_histTMajor); 
    normalize(_histTMinor); 
    normalize(_histOblateness); 

    normalize(_histSphericity); 
    normalize(_histAplanarity); 
    normalize(_histPlanarity); 

    normalize(_histHemiMassD); 
    normalize(_histHemiMassH); 
    normalize(_histHemiMassL); 

    normalize(_histHemiBroadW); 
    normalize(_histHemiBroadN); 
    normalize(_histHemiBroadT); 
    normalize(_histHemiBroadD); 

    normalize(_histCParam); 
    normalize(_histDParam); 

    normalize(_histDiffRate2Durham); 
    normalize(_histDiffRate2Jade); 
    normalize(_histDiffRate3Durham);
    normalize(_histDiffRate3Jade); 
    normalize(_histDiffRate4Durham);
    normalize(_histDiffRate4Jade); 
  }

}
