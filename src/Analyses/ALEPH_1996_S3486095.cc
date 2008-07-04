// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/ALEPH_1996_S3486095.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"

namespace Rivet {


  void ALEPH_1996_S3486095::analyze(const Event& e) {
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
    const ParticlePair& beams = applyProjection<Beam>(e, "Beams").getBeams();
    const double meanBeamMom = ( beams.first.getMomentum().vector3().mod() + 
                                 beams.second.getMomentum().vector3().mod() ) / 2.0;
    getLog() << Log::DEBUG << "Avg beam momentum = " << meanBeamMom << endl;

    // Thrusts
    getLog() << Log::DEBUG << "Calculating thrust" << endl;
    const Thrust& thrust = applyProjection<Thrust>(e, "Thrust");
    _hist1MinusT->fill(1 - thrust.thrust(), weight); 
    _histTMinor->fill(thrust.thrustMinor(), weight); 
    _histOblateness->fill(thrust.oblateness(), weight);

    // Jets
    #ifdef HAVE_JADE
    getLog() << Log::DEBUG << "Using FastJet JADE patch to make diff jet rate plots:" << endl;
    const FastJets& durjet = applyProjection<FastJets>(e, "DurhamJets");
    double y3 = durjet.getClusterSeq().exclusive_dmerge(2);
    _histY3->fill(-1. * std::log(y3), weight);

    #endif

    // Sphericities
    getLog() << Log::DEBUG << "Calculating sphericity" << endl;
    const Sphericity& sphericity = applyProjection<Sphericity>(e, "Sphericity");
    _histSphericity->fill(sphericity.sphericity(), weight); 
    _histAplanarity->fill(sphericity.aplanarity(), weight); 

    // C param
    getLog() << Log::DEBUG << "Calculating Parisi params" << endl;
    const ParisiTensor& parisi = applyProjection<ParisiTensor>(e, "Parisi");
    _histCParam->fill(parisi.C(), weight);

    // Hemispheres
    getLog() << Log::DEBUG << "Calculating hemisphere variables" << endl;
    const Hemispheres& hemi = applyProjection<Hemispheres>(e, "Hemispheres");
    _histHeavyJetMass->fill(hemi.getScaledM2high(), weight);

    // Iterate over all the charged final state particles.
    double Evis = 0.0;
    double rapt05 = 0.;
    double rapt10 = 0.;
    double rapt15 = 0.;
    double rapt20 = 0.;
    //int numChParticles = 0;
    getLog() << Log::DEBUG << "About to iterate over charged FS particles" << endl;
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      // Get momentum and energy of each particle.
      const Vector3 mom3 = p->getMomentum().vector3();
      const double energy = p->getMomentum().E();
      Evis += energy;
      _numChParticles += weight;

      // Scaled momenta.
      const double mom = mom3.mod();
      const double scaledMom = mom/meanBeamMom;
      const double logInvScaledMom = -std::log(scaledMom);
      _histLogScaledMom->fill(logInvScaledMom, weight); 
      _histScaledMom->fill(scaledMom, weight); 

      // Get momenta components w.r.t. thrust and sphericity.
      const double momT = dot(thrust.thrustAxis(), mom3);
      const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
      const double pToutS = dot(mom3, sphericity.sphericityMinorAxis());
      _histPtSIn->fill(fabs(pTinS/GeV), weight);
      _histPtSOut->fill(fabs(pToutS/GeV), weight);

      // Calculate rapidities w.r.t. thrust.
      const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
      _histRapidityT->fill(rapidityT, weight);
      if (std::fabs(rapidityT) <= .5)  {
          rapt05 += 1.;
      }
      if (std::fabs(rapidityT) <= 1.)  {
          rapt10 += 1.;
      }
      if (std::fabs(rapidityT) <= 1.5) {
          rapt15 += 1.;
      }
      if (std::fabs(rapidityT) <= 2.)  {
          rapt20 += 1.;
      } 

    }

    _histChMult->fill(numParticles, weight);

    _histMeanChMultRapt05->fill(_histMeanChMultRapt05->binMean(0), rapt05 * weight);
    _histMeanChMultRapt10->fill(_histMeanChMultRapt10->binMean(0), rapt10 * weight);
    _histMeanChMultRapt15->fill(_histMeanChMultRapt15->binMean(0), rapt15 * weight);
    _histMeanChMultRapt20->fill(_histMeanChMultRapt20->binMean(0), rapt20 * weight);
    _histMeanChMult->fill(_histMeanChMult->binMean(0), numParticles*weight);


    //// Final state of unstable particles to get particle spectra
    const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");
    for (ParticleVector::const_iterator p = ufs.particles().begin(); p != ufs.particles().end(); ++p) {
      const Vector3 mom3 = p->getMomentum().vector3();
      int id = abs(p->getPdgId());
      const double mom = mom3.mod();
      const double scaledMom = mom/meanBeamMom;
      switch (id) {
         case 22: 
            _histMultiPhoton->fill(-1.*std::log(scaledMom), weight);
            _weightedTotalNumPhoton += weight;
            break;
         case -321:
         case 321:
            _weightedTotalNumKPlus += weight;
            _histMultiKPlus->fill(scaledMom, weight);
            break;
         case 211:
         case -211:
            _histMultiPiPlus->fill(scaledMom, weight);
            _weightedTotalNumPiPlus += weight;
            break;
         case 2212:
         case -2212:
            _histMultiP->fill(scaledMom, weight);
            _weightedTotalNumP += weight;
            break;
         case 111:
            _histMultiPi0->fill(scaledMom, weight);
            _histMeanMultiPi0->fill(_histMeanMultiPi0->binMean(0), weight);
            _weightedTotalNumPi0 += weight;
            break;
         case 221:
            _histMultiEta->fill(scaledMom, weight);
            _histMeanMultiEta->fill(_histMeanMultiEta->binMean(0), weight);
            _weightedTotalNumEta += weight;
            break;
         case 331:
            _histMultiEtaPrime->fill(scaledMom, weight);
            _histMeanMultiEtaPrime->fill(_histMeanMultiEtaPrime->binMean(0), weight);
            _weightedTotalNumEtaPrime += weight;
            break;
         case 130: //klong
         case 310: //kshort
            _histMultiK0->fill(scaledMom, weight);
            _histMeanMultiK0->fill(_histMeanMultiK0->binMean(0), weight);
            _weightedTotalNumK0 += weight;
            break;
         case 113:
            _histMultiRho->fill(scaledMom, weight);
            _histMeanMultiRho->fill(_histMeanMultiRho->binMean(0), weight);
            _weightedTotalNumRho += weight;
            break;
         case 223:
            _histMeanMultiOmega->fill(_histMeanMultiOmega->binMean(0), weight);
            break;
         case 333:
            _histMultiPhi->fill(scaledMom, weight);
            _histMeanMultiPhi->fill(_histMeanMultiPhi->binMean(0), weight);
            _weightedTotalNumPhi += weight;
            break;
         case 313:
         case -313:
            _histMultiKStar892_0->fill(scaledMom, weight);
            _histMeanMultiKStar892_0->fill(_histMeanMultiKStar892_0->binMean(0), weight);
            _weightedTotalNumKStar892_0 += weight;
            break;
         case 323:
         case -323:
            _histMultiKStar892Plus->fill(scaledMom, weight);
            _histMeanMultiKStar892Plus->fill(_histMeanMultiKStar892Plus->binMean(0), weight);
            _weightedTotalNumKStar892Plus += weight;
            break;
         case 3122:
         case -3122:
            _histMultiLambda0->fill(scaledMom, weight);
            _histMeanMultiLambda0->fill(_histMeanMultiLambda0->binMean(0), weight);
            _weightedTotalNumLambda0 += weight;
            break;
         case 3212:
         case -3212:
            _histMeanMultiSigma0->fill(_histMeanMultiSigma0->binMean(0), weight);
         case 3312:
         case -3312:
            _histMeanMultiXiMinus->fill(scaledMom, weight);
            _histMeanMultiXiMinus->fill(_histMeanMultiXiMinus->binMean(0), weight);
            _weightedTotalNumXiMinus += weight;
            break;
         case 3114:
         case -3114: //maybe missing sigma(1385p13)
            _histMeanMultiSigma1385Plus->fill(_histMeanMultiSigma1385Plus->binMean(0), weight);
            _weightedTotalNumSigma1385Plus += weight;
            break;
         case 3324:
         case -3324:
            _histMeanMultiXi1530_0->fill(_histMeanMultiXi1530_0->binMean(0), weight);
            _weightedTotalNumXi1530_0 += weight;
            break;
         case 3334: 
            _histMultiOmegaMinus->fill(scaledMom, weight);
            _weightedTotalNumOmegaMinus += weight;
            break;
      }
    }

  }



  void ALEPH_1996_S3486095::init() {
    _histSphericity  = bookHistogram1D(1, 1, 1, "Sphericity, S (charged)");
    _histAplanarity  = bookHistogram1D(2, 1, 1, "Aplanarity, A (charged)");
    
    _hist1MinusT     = bookHistogram1D(3, 1, 1, "1-thrust, 1-T (charged)");
    _histTMinor      = bookHistogram1D(4, 1, 1, "Thrust minor, m (charged)");

    _histY3          = bookHistogram1D(5, 1, 1, "Two-jet resolution variable, Y3 (charged)");
    _histHeavyJetMass= bookHistogram1D(6, 1, 1, "Heavy Jet Mass (charged)"); 
    _histCParam      = bookHistogram1D(7, 1, 1, "C parameter (charged)");
    _histOblateness  = bookHistogram1D(8, 1, 1, "Oblateness = M - m (charged)");

    _histScaledMom   = bookHistogram1D(9, 1, 1, "Scaled momentum, x_p = |p|/|p_beam| (charged)");
    _histRapidityT   = bookHistogram1D(10, 1, 1, "Rapidity w.r.t. thrust axes, y_T (charged)");

    _histPtSIn       = bookHistogram1D(11, 1, 1, "In-plane p_T in GeV w.r.t. sphericity axes (charged)");
    _histPtSOut      = bookHistogram1D(12, 1, 1, "Out-of-plane p_T in GeV w.r.t. sphericity axes (charged)");

    _histLogScaledMom = bookHistogram1D(17, 1, 1, "Log of scaled momentum, log(1/x_p) (charged)");

    _histChMult       = bookHistogram1D(18, 1, 1, "Charged Multiplicity Distribution");
    _histMeanChMult   = bookHistogram1D(19, 1, 1, "Mean Charged Multiplicity");

    _histMeanChMultRapt05= bookHistogram1D(20, 1, 1, "Mean Charged Multiplicity for rapidity |Y| less than 0.5");
    _histMeanChMultRapt10= bookHistogram1D(21, 1, 1, "Mean Charged Multiplicity for rapidity |Y| less than 1.0");
    _histMeanChMultRapt15= bookHistogram1D(22, 1, 1, "Mean Charged Multiplicity for rapidity |Y| less than 1.5");
    _histMeanChMultRapt20= bookHistogram1D(23, 1, 1, "Mean Charged Multiplicity for rapidity |Y| less than 2.0");


    // particle spectra
    _histMultiPiPlus        = bookHistogram1D(25, 1, 1, "pi+/pi- multiplicity");
    _histMultiKPlus         = bookHistogram1D(26, 1, 1, "K+/K- multiplicity");
    _histMultiP             = bookHistogram1D(27, 1, 1, "p multiplicity");
    _histMultiPhoton        = bookHistogram1D(28, 1, 1, "photon multiplicity");
    _histMultiPi0           = bookHistogram1D(29, 1, 1, "pi0 multiplicity");
    _histMultiEta           = bookHistogram1D(30, 1, 1, "eta multiplicity");
    _histMultiEtaPrime      = bookHistogram1D(31, 1, 1, "etaprime multiplicity");
    _histMultiK0            = bookHistogram1D(32, 1, 1, "K0 multiplicity");
    _histMultiLambda0       = bookHistogram1D(33, 1, 1, "Lambda0 multiplicity");
    _histMultiXiMinus       = bookHistogram1D(34, 1, 1, "Xi- multiplicity");
    _histMultiSigma1385Plus = bookHistogram1D(35, 1, 1, "Sigma(1385)+/Sigma(1385)- multiplicity");
    _histMultiXi1530_0      = bookHistogram1D(36, 1, 1, "Xi(1530)0 multiplicity");
    _histMultiRho           = bookHistogram1D(37, 1, 1, "rho multiplicity");
    _histMultiOmegaMinus    = bookHistogram1D(38, 1, 1, "Omega- multiplicity");
    _histMultiKStar892_0    = bookHistogram1D(39, 1, 1, "K*(892)0 multiplicity");
    _histMultiPhi           = bookHistogram1D(40, 1, 1, "phi multiplicity");
    
    _histMultiKStar892Plus  = bookHistogram1D(43, 1, 1, "K*(892)+/K*(892)- multiplicity");
    
    // mean multiplicities 
    _histMeanMultiPi0           = bookHistogram1D(44, 1, 2, "Mean pi^0 multiplicity");
    _histMeanMultiEta           = bookHistogram1D(44, 1, 3, "Mean eta multiplicity");
    _histMeanMultiEtaPrime      = bookHistogram1D(44, 1, 4, "Mean etaprime multiplicity");
    _histMeanMultiK0            = bookHistogram1D(44, 1, 5, "Mean KS + KL multiplicity");
    _histMeanMultiRho           = bookHistogram1D(44, 1, 6, "Mean rho^0 multiplicity");
    _histMeanMultiOmega         = bookHistogram1D(44, 1, 7, "Mean omega^0 multiplicity");
    _histMeanMultiPhi           = bookHistogram1D(44, 1, 8, "Mean phi multiplicity");
    _histMeanMultiKStar892Plus  = bookHistogram1D(44, 1, 9, "Mean K*+ / K*- multiplicity");
    _histMeanMultiKStar892_0    = bookHistogram1D(44, 1, 10, "Mean K*0 / K*0bar multiplicity");
    _histMeanMultiLambda0       = bookHistogram1D(44, 1, 11, "Mean Lambda / Lambdabar multiplicity");
    _histMeanMultiSigma0        = bookHistogram1D(44, 1, 12, "Mean Sigma / Sigmabar multiplicity");
    _histMeanMultiXiMinus       = bookHistogram1D(44, 1, 13, "Mean Xi / Xibar multiplicity");
    _histMeanMultiSigma1385Plus = bookHistogram1D(44, 1, 14, "Mean Sigma(1385) / Sigma(1385)bar multiplicity");
    _histMeanMultiXi1530_0      = bookHistogram1D(44, 1, 15, "Mean Xi(1530) / Xi(1530)bar multiplicity");
  }



  // Finalize
  void ALEPH_1996_S3486095::finalize() { 
    // Normalize inclusive single particle distributions to the average number 
    // of charged particles per event.
    const double avgNumParts = _weightedTotalPartNum / sumOfWeights();

    normalize(_histPtSIn, avgNumParts);
    normalize(_histPtSOut, avgNumParts); 

    normalize(_histRapidityT, avgNumParts); 
    normalize(_histY3); 

    normalize(_histLogScaledMom, avgNumParts);
    normalize(_histScaledMom, avgNumParts); 

    // particle spectra
    scale(_histMultiPiPlus        ,1./sumOfWeights());
    scale(_histMultiKPlus         ,1./sumOfWeights());
    scale(_histMultiP             ,1./sumOfWeights());
    scale(_histMultiPhoton        ,1./sumOfWeights());
    scale(_histMultiPi0           ,1./sumOfWeights());
    scale(_histMultiEta           ,1./sumOfWeights());
    scale(_histMultiEtaPrime      ,1./sumOfWeights());
    scale(_histMultiK0            ,1./sumOfWeights());
    scale(_histMultiLambda0       ,1./sumOfWeights());
    scale(_histMultiXiMinus       ,1./sumOfWeights());
    scale(_histMultiSigma1385Plus ,1./sumOfWeights());
    scale(_histMultiXi1530_0      ,1./sumOfWeights());
    scale(_histMultiRho           ,1./sumOfWeights());
    scale(_histMultiOmegaMinus    ,1./sumOfWeights());
    scale(_histMultiKStar892_0    ,1./sumOfWeights());
    scale(_histMultiPhi           ,1./sumOfWeights());

    scale(_histMultiKStar892Plus  ,1./sumOfWeights());
    
    //normalize(_histMultiPiPlus        ,_weightedTotalNumPiPlus / sumOfWeights());
    //normalize(_histMultiKPlus         ,_weightedTotalNumKPlus/sumOfWeights());
    //normalize(_histMultiP             ,_weightedTotalNumP/sumOfWeights());
    //normalize(_histMultiPhoton            ,_weightedTotalNumPhoton/sumOfWeights());
    //normalize(_histMultiPi0           ,_weightedTotalNumPi0/sumOfWeights());
    //normalize(_histMultiEta           ,_weightedTotalNumEta/sumOfWeights());
    //normalize(_histMultiEtaPrime      ,_weightedTotalNumEtaPrime/sumOfWeights());
    //normalize(_histMultiK0            ,_weightedTotalNumK0/sumOfWeights());
    //normalize(_histMultiLambda0       ,_weightedTotalNumLambda0/sumOfWeights());
    //normalize(_histMultiXiMinus       ,_weightedTotalNumXiMinus/sumOfWeights());
    //normalize(_histMultiSigma1385Plus ,_weightedTotalNumSigma1385Plus/sumOfWeights());
    //normalize(_histMultiXi1530_0      ,_weightedTotalNumXi1530_0 /sumOfWeights());
    //normalize(_histMultiRho           ,_weightedTotalNumRho/sumOfWeights());
    //normalize(_histMultiOmegaMinus    ,_weightedTotalNumOmegaMinus/sumOfWeights());
    //normalize(_histMultiKStar892_0    ,_weightedTotalNumKStar892_0/sumOfWeights());
    //normalize(_histMultiPhi           ,_weightedTotalNumPhi/sumOfWeights());

    //normalize(_histMultiKStar892Plus  ,_weightedTotalNumKStar892Plus/sumOfWeights());

    // event shape
    normalize(_hist1MinusT); 
    normalize(_histTMinor); 
    normalize(_histOblateness); 

    normalize(_histSphericity); 
    normalize(_histAplanarity); 
    normalize(_histHeavyJetMass);  
    normalize(_histCParam); 
   

    // mean multiplicities 
    scale(_histChMult              , 2.0/sumOfWeights()); // taking into account the binwidth of 2 
    scale(_histMeanChMult          , 1.0/sumOfWeights());
    scale(_histMeanChMultRapt05    , 1.0/sumOfWeights());
    scale(_histMeanChMultRapt10    , 1.0/sumOfWeights());
    scale(_histMeanChMultRapt15    , 1.0/sumOfWeights());
    scale(_histMeanChMultRapt20    , 1.0/sumOfWeights());


    scale(_histMeanMultiPi0          , 1.0/sumOfWeights());
    scale(_histMeanMultiEta          , 1.0/sumOfWeights());
    scale(_histMeanMultiEtaPrime     , 1.0/sumOfWeights());
    scale(_histMeanMultiK0           , 1.0/sumOfWeights());
    scale(_histMeanMultiRho          , 1.0/sumOfWeights());
    scale(_histMeanMultiOmega        , 1.0/sumOfWeights());
    scale(_histMeanMultiPhi          , 1.0/sumOfWeights());
    scale(_histMeanMultiKStar892Plus , 1.0/sumOfWeights());
    scale(_histMeanMultiKStar892_0   , 1.0/sumOfWeights());
    scale(_histMeanMultiLambda0      , 1.0/sumOfWeights());
    scale(_histMeanMultiSigma0, 1.0/sumOfWeights());
    scale(_histMeanMultiXiMinus      , 1.0/sumOfWeights());
    scale(_histMeanMultiSigma1385Plus, 1.0/sumOfWeights());
    scale(_histMeanMultiXi1530_0     , 1.0/sumOfWeights());
  }

}
