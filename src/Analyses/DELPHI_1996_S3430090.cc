// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/DELPHI_1996_S3430090.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  DELPHI_1996_S3430090::DELPHI_1996_S3430090() 
  {
    setBeams(ELECTRON, POSITRON); 
    addProjection(Beam(), "Beams");
    const ChargedFinalState cfs;
    addProjection(cfs, "FS");
    addProjection(UnstableFinalState(), "UFS");
    addProjection(FastJets(cfs, FastJets::JADE, 0.7), "JadeJets");
    addProjection(FastJets(cfs, FastJets::DURHAM, 0.7), "DurhamJets");
    addProjection(Sphericity(cfs), "Sphericity");
    addProjection(ParisiTensor(cfs), "Parisi");
    const Thrust thrust(cfs);
    addProjection(thrust, "Thrust");
    addProjection(Hemispheres(thrust), "Hemispheres");
    _weightedTotalPartNum = 0;
    _passedCutWeightSum = 0;
  }



  void DELPHI_1996_S3430090::analyze(const Event& e) {
    // First, veto on leptonic events by requiring at least 4 charged FS particles
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const size_t numParticles = fs.particles().size();

    // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
    if (numParticles < 2) {
      getLog() << Log::DEBUG << "Failed leptonic event cut" << endl;
      vetoEvent;
    }
    getLog() << Log::DEBUG << "Passed leptonic event cut" << endl;
    const double weight = e.weight();
    _passedCutWeightSum += weight;
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

    foreach (const Particle& p, ufs.particles()) {
      int id = abs(p.pdgId());
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


  string DELPHI_1996_S3430090::dsigbyd(const string& x) {
    return "\\text{d}{\\sigma}/\\text{d}{" + x + "}";
  }
  string DELPHI_1996_S3430090::Ndsigbyd(const string& x) {
    return "N \\, " + dsigbyd(x);
  }
  string DELPHI_1996_S3430090::unitdsigbyd(const string& x) {
    return "N \\, " + dsigbyd(x);
  }
  string DELPHI_1996_S3430090::texmath(const string& foo) {
    return "$" + foo + "$";
  }


  void DELPHI_1996_S3430090::init() {
    _histPtTIn = 
      bookHistogram1D(1, 1, 1, "In-plane $p_\\perp$ in GeV w.r.t. thrust axes", 
                      "$p_\\perp^\\text{in}$ / GeV", texmath(Ndsigbyd("p_\\perp^\\text{in}")));
    _histPtTOut =
      bookHistogram1D(2, 1, 1, "Out-of-plane $p_\\perp$ in GeV w.r.t. thrust axes",
                      "$p_\\perp^\\text{out}$ / GeV", texmath(Ndsigbyd("p_\\perp^\\text{out}")));
    _histPtSIn =
      bookHistogram1D(3, 1, 1, "In-plane $p_\\perp$ in GeV w.r.t. sphericity axes",
                      "$p_\\perp^\\text{in}$ / GeV", texmath(Ndsigbyd("p_\\perp^\\text{in}")));
    _histPtSOut =
      bookHistogram1D(4, 1, 1, "Out-of-plane $p_\\perp$ in GeV w.r.t. sphericity axes",
                      "$p_\\perp^\\text{out}$ / GeV", texmath(Ndsigbyd("p_\\perp^\\text{out}")));
    
    _histRapidityT =
      bookHistogram1D(5, 1, 1, "Rapidity w.r.t. thrust axes, $y_T$",
                      "$y_T$", texmath(Ndsigbyd("y_T")));
    _histRapidityS =
      bookHistogram1D(6, 1, 1, "Rapidity w.r.t. sphericity axes, $y_S$",
                      "$y_S$", texmath(Ndsigbyd("y_S")));
    _histScaledMom =
      bookHistogram1D(7, 1, 1, "Scaled momentum, $x_p = |p|/|p_\\text{beam}|$",
                      "$x_p$", texmath(Ndsigbyd("x_p")));
    _histLogScaledMom =
      bookHistogram1D(8, 1, 1, "Log of scaled momentum, $\\log(1/x_p)$",
                      "$\\log(1/x_p)$", texmath(Ndsigbyd("\\log(1/x_p)")));
    
    _histPtTOutVsXp =
      bookProfile1D(9,  1, 1, "Mean out-of-plane $p_\\perp$ in GeV w.r.t. thrust axes vs. $x_p$",
                    "$x_p$", "$p_\\perp^\\text{out}$");
    _histPtVsXp =
      bookProfile1D(10, 1, 1, "Mean $p_\\perp$ in GeV vs. $x_p$",
                    "$x_p$", "$p_\\perp$");    

    _hist1MinusT =
      bookHistogram1D(11, 1, 1, "$1-\\text{Thrust}$", 
                      "$1-T$", texmath(unitdsigbyd("(1-T)")));
    _histTMajor =
      bookHistogram1D(12, 1, 1, "Thrust major, $M$",
                      "$M$", texmath(unitdsigbyd("M")));
    _histTMinor =
      bookHistogram1D(13, 1, 1, "Thrust minor, $m$",
                      "$m$", texmath(unitdsigbyd("m")));
    _histOblateness =
      bookHistogram1D(14, 1, 1, "Oblateness = $M - m$",
                      "$O$", texmath(unitdsigbyd("O")));
    
    _histSphericity =
      bookHistogram1D(15, 1, 1, "Sphericity, $S$",
                      "$S$", texmath(unitdsigbyd("S")));
    _histAplanarity =
      bookHistogram1D(16, 1, 1, "Aplanarity, $A$",
                      "$A$", texmath(unitdsigbyd("A")));
    _histPlanarity =
      bookHistogram1D(17, 1, 1, "Planarity, $P$",
                      "$P$", texmath(unitdsigbyd("P")));
    
    _histCParam =
      bookHistogram1D(18, 1, 1, "$C$ parameter",
                      "$C$", texmath(unitdsigbyd("C")));
    _histDParam =
      bookHistogram1D(19, 1, 1, "$D$ parameter",
                      "$D$", texmath(unitdsigbyd("D")));
    
    string y = "M_h^2/E_\\text{vis}^2";
    _histHemiMassH =
      bookHistogram1D(20, 1, 1, "Heavy hemisphere masses, $" + y + "$",
                      texmath(y), texmath(unitdsigbyd(y)));
    y = "M_l^2/E_\\text{vis}^2";
    _histHemiMassL =
      bookHistogram1D(21, 1, 1, "Light hemisphere masses, $" + y + "$",
                      texmath(y), texmath(unitdsigbyd(y)));
    y = "M_d^2/E_\\text{vis}^2";
    _histHemiMassD =
      bookHistogram1D(22, 1, 1, "Difference in hemisphere masses, $" + y + "$",
                      texmath(y), texmath(unitdsigbyd(y)));

    y = "B_\\text{max}";
    _histHemiBroadW =
      bookHistogram1D(23, 1, 1, "Wide hemisphere broadening, $" + y + "$",
                      texmath(y), texmath(unitdsigbyd(y)));
    y = "B_\\text{min}";
    _histHemiBroadN =
      bookHistogram1D(24, 1, 1, "Narrow hemisphere broadening, $B_\\text{min}$",
                      texmath(y), texmath(unitdsigbyd(y)));
    y = "B_\\text{sum}";
    _histHemiBroadT =
      bookHistogram1D(25, 1, 1, "Total hemisphere broadening, $B_\\text{sum}$",
                      texmath(y), texmath(unitdsigbyd(y)));
    y = "B_\\text{diff}";
    _histHemiBroadD =
      bookHistogram1D(26, 1, 1, "Difference in hemisphere broadening, $B_\\text{diff}$",
                      texmath(y), texmath(unitdsigbyd(y)));

    // Binned in y_cut
    y = "D_2^\\text{Durham}";
    _histDiffRate2Durham =
      bookHistogram1D(27, 1, 1, "Differential 2-jet rate with Durham algorithm, $" + y + "$",
                      texmath(y), texmath(unitdsigbyd(y)));                      
    y = "D_2^\\text{Jade}";
    _histDiffRate2Jade =
      bookHistogram1D(28, 1, 1, "Differential 2-jet rate with Jade algorithm, $D_2^\\text{Jade}$",
                      texmath(y), texmath(unitdsigbyd(y)));
    y = "D_3^\\text{Durham}";
    _histDiffRate3Durham =
      bookHistogram1D(29, 1, 1, "Differential 3-jet rate with Durham algorithm, $D_3^\\text{Durham}$",
                      texmath(y), texmath(unitdsigbyd(y)));
    y = "D_3^\\text{Jade}";
    _histDiffRate3Jade =
      bookHistogram1D(30, 1, 1, "Differential 3-jet rate with Jade algorithm, $D_3^\\text{Jade}$",
                      texmath(y), texmath(unitdsigbyd(y)));
    y = "D_4^\\text{Durham}";
    _histDiffRate4Durham =
      bookHistogram1D(31, 1, 1, "Differential 4-jet rate with Durham algorithm, $D_4^\\text{Durham}$",
                      texmath(y), texmath(unitdsigbyd(y)));
    y = "D_4^\\text{Jade}";
    _histDiffRate4Jade =
      bookHistogram1D(32, 1, 1, "Differential 4-jet rate with Jade algorithm, $D_4^\\text{Jade}$",
                      texmath(y), texmath(unitdsigbyd(y)));

    // Binned in cos(chi)
    string coschi = "$\\cos{\\chi}$";
    _histEEC =
      bookHistogram1D(33, 1, 1, "Energy-energy correlation, EEC", coschi, "EEC");
    _histAEEC =
      bookHistogram1D(34, 1, 1, "Asymmetry of the energy-energy correlation, AEEC", coschi, "AEEC");

    _histMultiCharged =
      bookHistogram1D(35, 1, 1, "Mean charged multiplicity", "", "Multiplicity");
    
    _histMultiPiPlus =
      bookHistogram1D(36, 1, 1, "Mean $\\pi^+/\\pi^-$ multiplicity", "", "Multiplicity");
    _histMultiPi0 =
      bookHistogram1D(36, 1, 2, "Mean $\\pi^0$ multiplicity", "", "Multiplicity");
    _histMultiKPlus =
      bookHistogram1D(36, 1, 3, "Mean $K^+/K^-$ multiplicity", "", "Multiplicity");
    _histMultiK0 =
      bookHistogram1D(36, 1, 4, "Mean $K^0$ multiplicity", "", "Multiplicity");
    _histMultiEta =
      bookHistogram1D(36, 1, 5, "Mean $\\eta$ multiplicity", "", "Multiplicity");
    _histMultiEtaPrime =
      bookHistogram1D(36, 1, 6, "Mean $\\eta'$ multiplicity", "", "Multiplicity");
    _histMultiDPlus =
      bookHistogram1D(36, 1, 7, "Mean $D^+$ multiplicity", "", "Multiplicity");
    _histMultiD0 =
      bookHistogram1D(36, 1, 8, "Mean $D^0$ multiplicity", "", "Multiplicity");
    _histMultiBPlus0 =
      bookHistogram1D(36, 1, 9, "Mean $B^+/B^-/B^0$ multiplicity", "", "Multiplicity");
    
    _histMultiF0 =
      bookHistogram1D(37, 1, 1, "Mean $f_0(980)$ multiplicity", "", "Multiplicity");
    
    _histMultiRho =
      bookHistogram1D(38, 1, 1, "Mean $\\rho$ multiplicity", "", "Multiplicity");
    _histMultiKStar892Plus =
      bookHistogram1D(38, 1, 2, "Mean $K^*(892)^+/K^*(892)^-$ multiplicity", "", "Multiplicity");
    _histMultiKStar892_0 =
      bookHistogram1D(38, 1, 3, "Mean $K^*(892)^0$ multiplicity", "", "Multiplicity");
    _histMultiPhi =
      bookHistogram1D(38, 1, 4, "Mean $\\phi$ multiplicity", "", "Multiplicity");
    _histMultiDStar2010Plus =
      bookHistogram1D(38, 1, 5, "Mean $D^*(2010)^+/D^*(2010)^-$ multiplicity", "", "Multiplicity");
    
    _histMultiF2 =
      bookHistogram1D(39, 1, 1, "Mean $f_2(1270)$ multiplicity", "", "Multiplicity");
    _histMultiK2Star1430_0 =
      bookHistogram1D(39, 1, 2, "Mean $K_2^*(1430)^0$ multiplicity", "", "Multiplicity");
    
    _histMultiP =
      bookHistogram1D(40, 1, 1, "Mean p multiplicity", "", "Multiplicity");
    _histMultiLambda0 =
      bookHistogram1D(40, 1, 2, "Mean $\\Lambda^0$ multiplicity", "", "Multiplicity");
    _histMultiXiMinus =
      bookHistogram1D(40, 1, 3, "Mean $\\Xi^-$ multiplicity", "", "Multiplicity");
    _histMultiOmegaMinus =
      bookHistogram1D(40, 1, 4, "Mean $\\Omega^-$ multiplicity", "", "Multiplicity");
    _histMultiDeltaPlusPlus =
      bookHistogram1D(40, 1, 5, "Mean $\\Delta(1232)^{++}$ multiplicity", "", "Multiplicity");
    _histMultiSigma1385Plus =
      bookHistogram1D(40, 1, 6, "Mean $\\Sigma(1385)^+/\\Sigma(1385)^-$ multiplicity", "", "Multiplicity");
    _histMultiXi1530_0 =
      bookHistogram1D(40, 1, 7, "Mean $\\Xi(1530)^0$ multiplicity", "", "Multiplicity");
    _histMultiLambdaB0 =
      bookHistogram1D(40, 1, 8, "Mean $\\Lambda_b^0$ multiplicity", "", "Multiplicity");
  }
  
  
  
  // Finalize
  void DELPHI_1996_S3430090::finalize() { 
    // Normalize inclusive single particle distributions to the average number 
    // of charged particles per event.
    const double avgNumParts = _weightedTotalPartNum / _passedCutWeightSum;

    normalize(_histPtTIn, avgNumParts);
    normalize(_histPtTOut, avgNumParts); 
    normalize(_histPtSIn, avgNumParts);
    normalize(_histPtSOut, avgNumParts); 

    normalize(_histRapidityT, avgNumParts); 
    normalize(_histRapidityS, avgNumParts); 

    normalize(_histLogScaledMom, avgNumParts);
    normalize(_histScaledMom, avgNumParts); 

    scale(_histEEC, 1.0/_passedCutWeightSum);
    scale(_histAEEC, 1.0/_passedCutWeightSum);
    scale(_histMultiCharged, 1.0/_passedCutWeightSum);

    scale(_histMultiPiPlus, 1.0/_passedCutWeightSum);
    scale(_histMultiPi0, 1.0/_passedCutWeightSum);
    scale(_histMultiKPlus, 1.0/_passedCutWeightSum);
    scale(_histMultiK0, 1.0/_passedCutWeightSum);
    scale(_histMultiEta, 1.0/_passedCutWeightSum);
    scale(_histMultiEtaPrime, 1.0/_passedCutWeightSum);
    scale(_histMultiDPlus, 1.0/_passedCutWeightSum);
    scale(_histMultiD0, 1.0/_passedCutWeightSum);
    scale(_histMultiBPlus0, 1.0/_passedCutWeightSum);

    scale(_histMultiF0, 1.0/_passedCutWeightSum);

    scale(_histMultiRho, 1.0/_passedCutWeightSum);
    scale(_histMultiKStar892Plus, 1.0/_passedCutWeightSum);
    scale(_histMultiKStar892_0, 1.0/_passedCutWeightSum);
    scale(_histMultiPhi, 1.0/_passedCutWeightSum);
    scale(_histMultiDStar2010Plus, 1.0/_passedCutWeightSum);

    scale(_histMultiF2, 1.0/_passedCutWeightSum);
    scale(_histMultiK2Star1430_0, 1.0/_passedCutWeightSum);

    scale(_histMultiP, 1.0/_passedCutWeightSum);
    scale(_histMultiLambda0, 1.0/_passedCutWeightSum);
    scale(_histMultiXiMinus, 1.0/_passedCutWeightSum);
    scale(_histMultiOmegaMinus, 1.0/_passedCutWeightSum);
    scale(_histMultiDeltaPlusPlus, 1.0/_passedCutWeightSum);
    scale(_histMultiSigma1385Plus, 1.0/_passedCutWeightSum);
    scale(_histMultiXi1530_0, 1.0/_passedCutWeightSum);
    scale(_histMultiLambdaB0, 1.0/_passedCutWeightSum);

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
