// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/DELPHI_1996_S3430090.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"


namespace Rivet {


  void DELPHI_1996_S3430090::analyze(const Event& e) {
    Log& log = getLog();

    // First, veto on leptonic events by examining the signal vertex
    // NB. Doesn't work --- signal vertex not reliable.
    // GenVertex* sigVtx = e.genEvent().signal_process_vertex();
    // size_t count = 0;
    // for (HepMC::GenVertex::particles_out_const_iterator sigp = sigVtx->particles_out_const_begin();
    //      sigp != sigVtx->particles_out_const_end(); ++sigp) {
    //   if (*sigp) {
    //     ++count;
    //     log << Log::INFO << "Sig part ID = " << (*sigp)->pdg_id() << endl;
    //   }
    // }
    // log << Log::INFO << "Num sig particles = " << count << endl;

    // First, veto on leptonic events by requiring at least 4 charged FS particles
    /// @todo Work out why this creates a segfault if it's run before the beams projection...
    const FinalState& fs = e.applyProjection(_cnfsproj);

    // FIXME: This lepton veto cut is merely a hack, but not a solution.
    // We should really tell the generator to produces only hadronic events.
    // Cutting on numCharged>=4 is not sufficient, and I have no idea how
    // good or bad numCharged>=6 really is. And the way we fix the normalisation
    // of single particle distributions is not nice, either.
    size_t numCharged = 0;
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      if (PID::threeCharge(p->getPdgId())) ++numCharged;
    }
    const string result = (numCharged < 4) ? "Failed" : "Passed";
    log << Log::DEBUG << result << " leptonic event cut" << endl;
    if (numCharged < 6) {
      _sumOfRejectedWeights+=e.weight();
      return;
    }

    // Get beams and average beam momentum
    const ParticlePair& beams = e.applyProjection(_beamsproj).getBeams();
    const double meanBeamMom = ( beams.first.getMomentum().vector3().mod() + 
                                 beams.second.getMomentum().vector3().mod() ) / 2.0;
    log << Log::DEBUG << "Avg beam momentum = " << meanBeamMom << endl;

    // Get event weight for histo filling
    const double weight = e.weight();

    // Thrusts
    log << Log::DEBUG << "Calculating thrust" << endl;
    const Thrust& thrustC = e.applyProjection(_cthrustproj);
    _hist1MinusTC->fill(1 - thrustC.thrust(), weight); 
    _hist1MinusTC->fill(1 - thrustC.thrust(), weight); 
    _histTMajorC->fill(thrustC.thrustMajor(), weight); 
    _histTMinorC->fill(thrustC.thrustMinor(), weight); 
    _histOblatenessC->fill(thrustC.oblateness(), weight);

    log << Log::DEBUG << "Calculating CN thrust" << endl;
    const Thrust& thrustCN = e.applyProjection(_cnthrustproj);
    _hist1MinusTCN->fill(1 - thrustCN.thrust(), weight); 
    _histTMajorCN->fill(thrustCN.thrustMajor(), weight); 
    _histTMinorCN->fill(thrustCN.thrustMinor(), weight); 
    _histOblatenessCN->fill(thrustCN.oblateness(), weight);

    // Jets
#ifdef __HAVE_JADE
     const FastJets& durjetC = e.applyProjection(_cdurjetproj);
     _histDiffRate2DurhamC->fill(durjetC.getClusterSeq().exclusive_dmerge(2), weight); 
     _histDiffRate3DurhamC->fill(durjetC.getClusterSeq().exclusive_dmerge(3), weight); 
     _histDiffRate4DurhamC->fill(durjetC.getClusterSeq().exclusive_dmerge(4), weight); 
     const FastJets& durjetCN = e.applyProjection(_cndurjetproj);
     _histDiffRate2DurhamCN->fill(durjetCN.getClusterSeq().exclusive_dmerge(2), weight); 
     _histDiffRate3DurhamCN->fill(durjetCN.getClusterSeq().exclusive_dmerge(3), weight); 
     _histDiffRate4DurhamCN->fill(durjetCN.getClusterSeq().exclusive_dmerge(4), weight); 
     const FastJets& jadejetC = e.applyProjection(_cjadejetproj);
     _histDiffRate2JadeC->fill(jadejetC.getClusterSeq().exclusive_dmerge(2), weight); 
     _histDiffRate3JadeC->fill(jadejetC.getClusterSeq().exclusive_dmerge(3), weight); 
     _histDiffRate4JadeC->fill(jadejetC.getClusterSeq().exclusive_dmerge(4), weight); 
     const FastJets& jadejetCN = e.applyProjection(_cnjadejetproj);
     _histDiffRate2JadeCN->fill(jadejetCN.getClusterSeq().exclusive_dmerge(2), weight); 
     _histDiffRate3JadeCN->fill(jadejetCN.getClusterSeq().exclusive_dmerge(3), weight); 
     _histDiffRate4JadeCN->fill(jadejetCN.getClusterSeq().exclusive_dmerge(4), weight); 
#endif

    // Sphericities
    log << Log::DEBUG << "Calculating sphericity" << endl;
    const Sphericity& sphericityC = e.applyProjection(_cspherproj);
    _histSphericityC->fill(sphericityC.sphericity(), weight); 
    _histAplanarityC->fill(sphericityC.aplanarity(), weight); 
    _histPlanarityC->fill(sphericityC.planarity(), weight); 

    const Sphericity& sphericityCN = e.applyProjection(_cnspherproj);
    _histSphericityCN->fill(sphericityCN.sphericity(), weight); 
    _histAplanarityCN->fill(sphericityCN.aplanarity(), weight); 

    // C & D params
    log << Log::DEBUG << "Calculating Parisi params" << endl;
    const ParisiTensor& parisiC = e.applyProjection(_cparisiproj);    
    _histCParamC->fill(parisiC.C(), weight); 
    _histDParamC->fill(parisiC.D(), weight); 

    const ParisiTensor& parisiCN = e.applyProjection(_cnparisiproj);    
    _histCParamCN->fill(parisiCN.C(), weight); 
    _histDParamCN->fill(parisiCN.D(), weight); 

    // Hemispheres
    log << Log::DEBUG << "Calculating hemisphere variables" << endl;
    const Hemispheres& hemiC = e.applyProjection(_chemiproj);
    _histHemiMassHC->fill(hemiC.getScaledM2high(), weight); 
    _histHemiMassLC->fill(hemiC.getScaledM2low(), weight); 
    _histHemiMassDC->fill(hemiC.getScaledM2diff(), weight); 
    _histHemiBroadWC->fill(hemiC.getBmax(), weight); 
    _histHemiBroadNC->fill(hemiC.getBmin(), weight); 
    _histHemiBroadTC->fill(hemiC.getBsum(), weight); 
    _histHemiBroadDC->fill(hemiC.getBdiff(), weight); 

    const Hemispheres& hemiCN = e.applyProjection(_cnhemiproj);
    _histHemiMassHCN->fill(hemiCN.getScaledM2high(), weight); 
    _histHemiMassLCN->fill(hemiCN.getScaledM2low(), weight); 
    _histHemiMassDCN->fill(hemiCN.getScaledM2diff(), weight); 
    _histHemiBroadWCN->fill(hemiCN.getBmax(), weight); 
    _histHemiBroadNCN->fill(hemiCN.getBmin(), weight); 
    _histHemiBroadTCN->fill(hemiCN.getBsum(), weight); 
    _histHemiBroadDCN->fill(hemiCN.getBdiff(), weight); 


    // Iterate over all the charged final state particles.
    const FinalState& fsC = e.applyProjection(_cfsproj);
    double Evis = 0.;
    double Evis2 = 0.;
    log << Log::DEBUG << "About to iterate over charged FS particles" << endl;
    for (ParticleVector::const_iterator p = fsC.particles().begin(); p != fsC.particles().end(); ++p) {
      // Get momentum and energy of each particle.
      const Vector3 mom3 = p->getMomentum().vector3();
      const double energy = p->getMomentum().E();
      Evis += energy;

      // Scaled momenta.
      const double mom = mom3.mod();
      const double scaledMom = mom/meanBeamMom;
      const double logInvScaledMom = -std::log(scaledMom);
      _histLogScaledMom->fill(logInvScaledMom, weight); 
      _histScaledMom->fill(scaledMom, weight); 

      // Get momenta components w.r.t. thrust and sphericity.
      const double momT = dot(thrustC.thrustAxis(), mom3);
      const double momS = dot(sphericityC.sphericityAxis(), mom3);
      const double pTinT = dot(mom3, thrustC.thrustMajorAxis());
      const double pToutT = dot(mom3, thrustC.thrustMinorAxis());
      const double pTinS = dot(mom3, sphericityC.sphericityMajorAxis());
      const double pToutS = dot(mom3, sphericityC.sphericityMinorAxis());
      const double pT = sqrt(pow(pTinT, 2) + pow(pToutT, 2));
      _histPtTInC->fill(fabs(pTinT/GeV), weight);
      _histPtTOutC->fill(fabs(pToutT/GeV), weight);
      _histPtSInC->fill(fabs(pTinS/GeV), weight);
      _histPtSOutC->fill(fabs(pToutS/GeV), weight);
      _histPtVsXp->fill(scaledMom, fabs(pT/GeV), weight);
      _histPtTOutVsXp->fill(scaledMom, fabs(pToutT/GeV), weight);

      // Calculate rapidities w.r.t. thrust and sphericity.
      const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
      const double rapidityS = 0.5 * std::log((energy + momS) / (energy - momS));
      _histRapidityTC->fill(rapidityT, weight); 
      _histRapiditySC->fill(rapidityS, weight); 
    }
    Evis2 = Evis*Evis;
    for (ParticleVector::const_iterator p_i = fsC.particles().begin(); p_i != fsC.particles().end(); ++p_i) {
      for (ParticleVector::const_iterator p_j = p_i; p_j != fsC.particles().end(); ++p_j) {
        if (p_i == p_j) continue;
        const Vector3 mom3_i = p_i->getMomentum().vector3();
        const Vector3 mom3_j = p_j->getMomentum().vector3();
        const double energy_i = p_i->getMomentum().E();
        const double energy_j = p_j->getMomentum().E();
        const double cosij = dot(mom3_i.unit(), mom3_j.unit());
        const double eec = (energy_i*energy_j) / Evis2;
        _histEEC->fill(cosij, eec*weight);
        _histAEEC->fill( cosij,  eec*weight);
        _histAEEC->fill(-cosij, -eec*weight);
      }
    }
    const size_t numParticlesC = fsC.particles().size();
    _weightedTotalPartNumC += numParticlesC * weight;


    // Iterate over all the charged and neutral final state particles.
    const FinalState& fsCN = e.applyProjection(_cnfsproj);
    log << Log::DEBUG << "About to iterate over charged+neutral FS particles" << endl;
    for (ParticleVector::const_iterator p = fsCN.particles().begin(); p != fsCN.particles().end(); ++p) {
      // Get momentum and energy of each particle.
      const Vector3 mom3 = p->getMomentum().vector3();
      const double energy = p->getMomentum().E();

      // Get momenta components w.r.t. thrust and sphericity.
      const double momT = dot(thrustCN.thrustAxis(), mom3);
      const double momS = dot(sphericityCN.sphericityAxis(), mom3);
      const double pTinT = dot(mom3, thrustCN.thrustMajorAxis());
      const double pToutT = dot(mom3, thrustCN.thrustMinorAxis());
      const double pTinS = dot(mom3, sphericityCN.sphericityMajorAxis());
      const double pToutS = dot(mom3, sphericityCN.sphericityMinorAxis());
      _histPtTInCN->fill(fabs(pTinT/GeV), weight);
      _histPtTOutCN->fill(fabs(pToutT/GeV), weight);
      _histPtSInCN->fill(fabs(pTinS/GeV), weight);
      _histPtSOutCN->fill(fabs(pToutS/GeV), weight);

      // Calculate rapidities w.r.t. thrust and sphericity.
      const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
      const double rapidityS = 0.5 * std::log((energy + momS) / (energy - momS));
      _histRapidityTCN->fill(rapidityT, weight); 
      _histRapiditySCN->fill(rapidityS, weight); 
    }
    const size_t numParticlesCN = fsCN.particles().size();
    _weightedTotalPartNumCN += numParticlesCN * weight;
  }



  void DELPHI_1996_S3430090::init() {
    _histPtTInC       = bookHistogram1D(1, 1, 1, "In-plane p_T in GeV w.r.t. thrust axes (charged)");
    _histPtTInCN      = bookHistogram1D(1, 1, 2, "In-plane p_T in GeV w.r.t. thrust axes (charged and neutral)");
    _histPtTOutC      = bookHistogram1D(2, 1, 1, "Out-of-plane p_T in GeV w.r.t. thrust axes (charged)");
    _histPtTOutCN     = bookHistogram1D(2, 1, 2, "Out-of-plane p_T in GeV w.r.t. thrust axes (charged and neutral)");
    _histPtSInC       = bookHistogram1D(3, 1, 1, "In-plane p_T in GeV w.r.t. sphericity axes (charged)");
    _histPtSInCN      = bookHistogram1D(3, 1, 2, "In-plane p_T in GeV w.r.t. sphericity axes (charged and neutral)");
    _histPtSOutC      = bookHistogram1D(4, 1, 1, "Out-of-plane p_T in GeV w.r.t. sphericity axes (charged)");
    _histPtSOutCN     = bookHistogram1D(4, 1, 2, "Out-of-plane p_T in GeV w.r.t. sphericity axes (charged and neutral)");

    _histRapidityTC   = bookHistogram1D(5, 1, 1, "Rapidity w.r.t. thrust axes, y_T (charged)");
    _histRapidityTCN  = bookHistogram1D(5, 1, 2, "Rapidity w.r.t. thrust axes, y_T (charged and neutral)");
    _histRapiditySC   = bookHistogram1D(6, 1, 1, "Rapidity w.r.t. sphericity axes, y_S (charged)");
    _histRapiditySCN  = bookHistogram1D(6, 1, 2, "Rapidity w.r.t. sphericity axes, y_S (charged and neutral)");

    _histScaledMom    = bookHistogram1D(7, 1, 1, "Scaled momentum, x_p = |p|/|p_beam| (charged)");
    _histLogScaledMom = bookHistogram1D(8, 1, 1, "Log of scaled momentum, log(1/x_p) (charged)");

    _histPtTOutVsXp   = bookProfile1D(9,  1, 1, "Mean out-of-plane p_T in GeV w.r.t. thrust axes vs. x_p (charged)"); // binned in Xp
    _histPtVsXp       = bookProfile1D(10, 1, 1, "Mean p_T in GeV vs. x_p (charged)"); // binned in Xp

    _hist1MinusTC     = bookHistogram1D(11, 1, 1, "1-thrust, 1-T (charged)");
    _hist1MinusTCN    = bookHistogram1D(11, 1, 2, "1-thrust, 1-T (charged and neutral)");
    _histTMajorC      = bookHistogram1D(12, 1, 1, "Thrust major, M (charged)");
    _histTMajorCN     = bookHistogram1D(12, 1, 2, "Thrust major, M (charged and neutral)");
    _histTMinorC      = bookHistogram1D(13, 1, 1, "Thrust minor, m (charged)");
    _histTMinorCN     = bookHistogram1D(13, 1, 2, "Thrust minor, m (charged and neutral)");
    _histOblatenessC  = bookHistogram1D(14, 1, 1, "Oblateness = M - m (charged)");
    _histOblatenessCN = bookHistogram1D(14, 1, 2, "Oblateness = M - m (charged and neutral)");

    _histSphericityC  = bookHistogram1D(15, 1, 1, "Sphericity, S (charged)");
    _histSphericityCN = bookHistogram1D(15, 1, 2, "Sphericity, S (charged and neutral)");
    _histAplanarityC  = bookHistogram1D(16, 1, 1, "Aplanarity, A (charged)");
    _histAplanarityCN = bookHistogram1D(16, 1, 2, "Aplanarity, A (charged and neutral)");
    _histPlanarityC   = bookHistogram1D(17, 1, 1, "Planarity, P (charged)");
    _histPlanarityCN  = bookHistogram1D(17, 1, 2, "Planarity, P (charged and neutral)");

    _histCParamC      = bookHistogram1D(18, 1, 1, "C parameter (charged)");
    _histCParamCN     = bookHistogram1D(18, 1, 2, "C parameter (charged and neutral)");
    _histDParamC      = bookHistogram1D(19, 1, 1, "D parameter (charged)");
    _histDParamCN     = bookHistogram1D(19, 1, 2, "D parameter (charged and neutral)");

    _histHemiMassHC   = bookHistogram1D(20, 1, 1, "Heavy hemisphere masses, M_h^2/E_vis^2 (charged)");
    _histHemiMassHCN  = bookHistogram1D(20, 1, 2, "Heavy hemisphere masses, M_h^2/E_vis^2 (charged and neutral)");
    _histHemiMassLC   = bookHistogram1D(21, 1, 1, "Light hemisphere masses, M_l^2/E_vis^2 (charged)");
    _histHemiMassLCN  = bookHistogram1D(21, 1, 2, "Light hemisphere masses, M_l^2/E_vis^2 (charged and neutral)");
    _histHemiMassDC   = bookHistogram1D(22, 1, 1, "Difference in hemisphere masses, M_d^2/E_vis^2 (charged)");
    _histHemiMassDCN  = bookHistogram1D(22, 1, 2, "Difference in hemisphere masses, M_d^2/E_vis^2 (charged and neutral)");

    _histHemiBroadWC  = bookHistogram1D(23, 1, 1, "Wide hemisphere broadening, B_max (charged)");
    _histHemiBroadWCN = bookHistogram1D(23, 1, 2, "Wide hemisphere broadening, B_max (charged and neutral)");
    _histHemiBroadNC  = bookHistogram1D(24, 1, 1, "Narrow hemisphere broadening, B_min (charged)");
    _histHemiBroadNCN = bookHistogram1D(24, 1, 2, "Narrow hemisphere broadening, B_min (charged and neutral)");
    _histHemiBroadTC  = bookHistogram1D(25, 1, 1, "Total hemisphere broadening, B_sum (charged)");
    _histHemiBroadTCN = bookHistogram1D(25, 1, 2, "Total hemisphere broadening, B_sum (charged and neutral)");
    _histHemiBroadDC  = bookHistogram1D(26, 1, 1, "Difference in hemisphere broadening, B_diff (charged)");
    _histHemiBroadDCN = bookHistogram1D(26, 1, 2, "Difference in hemisphere broadening, B_diff (charged and neutral)");

    _histDiffRate2DurhamC  = bookHistogram1D(27, 1, 1, "Differential 2-jet rate with Durham algorithm, D_2^Durham (charged)"); // binned in y_cut
    _histDiffRate2DurhamCN = bookHistogram1D(27, 1, 2, "Differential 2-jet rate with Durham algorithm, D_2^Durham (charged and neutral)"); // binned in y_cut
    _histDiffRate2JadeC    = bookHistogram1D(28, 1, 1, "Differential 2-jet rate with Jade algorithm, D_2^Jade (charged)"); // binned in y_cut
    _histDiffRate2JadeCN   = bookHistogram1D(28, 1, 2, "Differential 2-jet rate with Jade algorithm, D_2^Jade (charged and neutral)"); // binned in y_cut
    _histDiffRate3DurhamC  = bookHistogram1D(29, 1, 1, "Differential 3-jet rate with Durham algorithm, D_3^Durham (charged)"); // binned in y_cut
    _histDiffRate3DurhamCN = bookHistogram1D(29, 1, 2, "Differential 3-jet rate with Durham algorithm, D_3^Durham (charged and neutral)"); // binned in y_cut
    _histDiffRate3JadeC    = bookHistogram1D(30, 1, 1, "Differential 3-jet rate with Jade algorithm, D_3^Jade (charged)"); // binned in y_cut
    _histDiffRate3JadeCN   = bookHistogram1D(30, 1, 2, "Differential 3-jet rate with Jade algorithm, D_3^Jade (charged and neutral)"); // binned in y_cut
    _histDiffRate4DurhamC  = bookHistogram1D(31, 1, 1, "Differential 4-jet rate with Durham algorithm, D_4^Durham (charged)"); // binned in y_cut
    _histDiffRate4DurhamCN = bookHistogram1D(31, 1, 2, "Differential 4-jet rate with Durham algorithm, D_4^Durham (charged and neutral)"); // binned in y_cut
    _histDiffRate4JadeC    = bookHistogram1D(32, 1, 1, "Differential 4-jet rate with Jade algorithm, D_4^Jade (charged)"); // binned in y_cut
    _histDiffRate4JadeCN   = bookHistogram1D(32, 1, 2, "Differential 4-jet rate with Jade algorithm, D_4^Jade (charged and neutral)"); // binned in y_cut

    _histEEC  = bookHistogram1D(33, 1, 1, "Energy-energy correlation, EEC (charged)"); // binned in cos(chi)
    _histAEEC = bookHistogram1D(34, 1, 1, "Asymmetry of the energy-energy correlation, AEEC (charged)"); // binned in cos(chi)


    // Identified particle distributions
    // ?
  }



  // Finalize
  void DELPHI_1996_S3430090::finalize() { 
    // Normalize inclusive single particle distributions to the average number 
    // of (charged / charged+neutral) particles per event.
    const double sumOfAcceptedWeights = sumOfWeights() - _sumOfRejectedWeights;
    const double avgNumPartsC = _weightedTotalPartNumC / sumOfAcceptedWeights;
    const double avgNumPartsCN = _weightedTotalPartNumCN / sumOfAcceptedWeights;
    
    normalize(_histPtTInC, avgNumPartsC);
    normalize(_histPtTInCN, avgNumPartsCN);
    normalize(_histPtTOutC, avgNumPartsC); 
    normalize(_histPtTOutCN, avgNumPartsCN); 
    normalize(_histPtSInC, avgNumPartsC);
    normalize(_histPtSInCN, avgNumPartsCN);
    normalize(_histPtSOutC, avgNumPartsC); 
    normalize(_histPtSOutCN, avgNumPartsCN); 
    
    normalize(_histRapidityTC, avgNumPartsC); 
    normalize(_histRapidityTCN, avgNumPartsCN);
    normalize(_histRapiditySC, avgNumPartsC); 
    normalize(_histRapiditySCN, avgNumPartsCN);
    
    normalize(_histLogScaledMom, avgNumPartsC);
    normalize(_histScaledMom, avgNumPartsC); 

    scale(_histEEC, 1.0/sumOfAcceptedWeights);
    scale(_histAEEC, 1.0/sumOfAcceptedWeights);

    normalize(_hist1MinusTC); 
    normalize(_hist1MinusTCN); 
    normalize(_histTMajorC); 
    normalize(_histTMajorCN); 
    normalize(_histTMinorC); 
    normalize(_histTMinorCN); 
    normalize(_histOblatenessC); 
    normalize(_histOblatenessCN); 
    
    normalize(_histSphericityC); 
    normalize(_histSphericityCN); 
    normalize(_histAplanarityC); 
    normalize(_histAplanarityCN); 
    normalize(_histPlanarityC); 
    normalize(_histPlanarityCN); 
    
    normalize(_histHemiMassDC); 
    normalize(_histHemiMassDCN); 
    normalize(_histHemiMassHC); 
    normalize(_histHemiMassHCN); 
    normalize(_histHemiMassLC); 
    normalize(_histHemiMassLCN); 
    
    normalize(_histHemiBroadWC); 
    normalize(_histHemiBroadWCN); 
    normalize(_histHemiBroadNC); 
    normalize(_histHemiBroadNCN); 
    normalize(_histHemiBroadTC); 
    normalize(_histHemiBroadTCN); 
    normalize(_histHemiBroadDC); 
    normalize(_histHemiBroadDCN); 
    
    normalize(_histCParamC); 
    normalize(_histCParamCN); 
    normalize(_histDParamC); 
    normalize(_histDParamCN); 
    
    normalize(_histDiffRate2DurhamC); 
    normalize(_histDiffRate2DurhamCN); 
    normalize(_histDiffRate2JadeC); 
    normalize(_histDiffRate2JadeCN);
    normalize(_histDiffRate3DurhamC);
    normalize(_histDiffRate3DurhamCN);
    normalize(_histDiffRate3JadeC); 
    normalize(_histDiffRate3JadeCN);
    normalize(_histDiffRate4DurhamC);
    normalize(_histDiffRate4DurhamCN);
    normalize(_histDiffRate4JadeC); 
    normalize(_histDiffRate4JadeCN); 
  }


}
