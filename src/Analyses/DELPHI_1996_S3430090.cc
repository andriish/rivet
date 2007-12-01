// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/DELPHI_1996_S3430090.hh"
#include "HepPDT/ParticleID.hh"

using namespace AIDA;
using namespace HepMC;


namespace Rivet {


  void DELPHI_1996_S3430090::analyze(const Event& e) {
    Log& log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;

    // Get beams and average beam momentum
    const ParticlePair& beams = e.applyProjection(_beamsproj).getBeams();
    const double meanBeamMom = ( beams.first.getMomentum().vector3().mod() + 
                                 beams.second.getMomentum().vector3().mod() ) / 2.0;

    // Get event weight for histo filling
    const double weight = e.weight();

    // Thrusts
    const Thrust& thrustC = e.applyProjection(_cthrustproj);
    _hist1MinusTC->fill(1 - thrustC.thrust(), weight); 
    _histTMajorC->fill(thrustC.thrustMajor(), weight); 
    _histTMinorC->fill(thrustC.thrustMinor(), weight); 
    _histOblatenessC->fill(thrustC.oblateness(), weight);

    const Thrust& thrustCN = e.applyProjection(_cnthrustproj);
    _hist1MinusTCN->fill(1 - thrustCN.thrust(), weight); 
    _histTMajorCN->fill(thrustCN.thrustMajor(), weight); 
    _histTMinorCN->fill(thrustCN.thrustMinor(), weight); 
    _histOblatenessCN->fill(thrustCN.oblateness(), weight);

    // Jets
    //const DurhamJet& durjetC = e.applyProjection(_cdurproj);
    //_histDiffRate2Durham->fill(, weight); 
    //_histDiffRate3Durham->fill(, weight);
    //_histDiffRate4Durham->fill(, weight);
    //const DurhamJet& durjetCN = e.applyProjection(_cndurproj);
    //_histDiffRate2DurhamCN->fill(, weight); 
    //_histDiffRate3DurhamCN->fill(, weight);
    //_histDiffRate4DurhamCN->fill(, weight);
    //const JadeJet& durjetC = e.applyProjection(_cjadeproj);
    //_histDiffRate2Jade->fill(, weight); 
    //_histDiffRate3Jade->fill(, weight); 
    //_histDiffRate4Jade->fill(, weight);
    //const JadeJet& durjetCN = e.applyProjection(_cnjadeproj);
    //_histDiffRate2JadeCN->fill(, weight); 
    //_histDiffRate3JadeCN->fill(, weight);
    //_histDiffRate4JadeCN->fill(, weight);

    // Sphericities
    const Sphericity& sphericityC = e.applyProjection(_cspherproj);
    _histSphericityC->fill(sphericityC.sphericity(), weight); 
    _histAplanarityC->fill(sphericityC.aplanarity(), weight); 
    _histPlanarityC->fill(sphericityC.planarity(), weight); 

    const Sphericity& sphericityCN = e.applyProjection(_cnspherproj);
    _histSphericityCN->fill(sphericityCN.sphericity(), weight); 
    _histAplanarityCN->fill(sphericityCN.aplanarity(), weight); 

    // C & D params
    const ParisiTensor& parisiC = e.applyProjection(_cparisiproj);    
    _histCParamC->fill(parisiC.C(), weight); 
    _histDParamC->fill(parisiC.D(), weight); 

    const ParisiTensor& parisiCN = e.applyProjection(_cnparisiproj);    
    _histCParamCN->fill(parisiCN.C(), weight); 
    _histDParamCN->fill(parisiCN.D(), weight); 

    // Hemispheres
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

    //const EEC& eec = e.applyProjection(_eecproj);
    //_histEEC->fill(, weight); 
    //_histAEEC->fill(, weight); 


    // Iterate over all the final state particles.
    double totalPtInT(0), totalPtOutT(0);
    const FinalState& cfsC = e.applyProjection(_cfsproj);
    const size_t numParticles = cfsC.particles().size();
    for (ParticleVector::const_iterator p = cfsC.particles().begin(); p != cfsC.particles().end(); ++p) {
      /// @todo Add missing CN variants

      // Get momentum and energy of each particle.
      const Vector3 mom3 = p->getMomentum().vector3();
      const double mom = mom3.mod();
      const double energy = p->getMomentum().E();

      // Calculate scaled momenta.
      const double scaledMom = mom/meanBeamMom;
      const double logInvScaledMom = -log10(scaledMom);

      // Get momenta components w.r.t. thrust and sphericity.
      const double momT = dot(thrustC.thrustAxis(), mom3);
      const double momS = dot(sphericityC.sphericityAxis(), mom3);
      const double pTinT = dot(mom3, thrustC.thrustMajorAxis());
      const double pToutT = dot(mom3, thrustC.thrustMinorAxis());
      const double pTinS = dot(mom3, sphericityC.sphericityMinorAxis());
      const double pToutS = dot(mom3, sphericityC.sphericityMinorAxis());
      totalPtInT += pTinT; 
      totalPtOutT += pToutT;

      // Calculate rapidities w.r.t. thrust and sphericity.
      const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
      const double rapidityS = 0.5 * std::log((energy + momS) / (energy - momS));

      // Fill histograms.
      _histLogScaledMom->fill(logInvScaledMom, weight); 
      _histScaledMom->fill(scaledMom, weight); 
      _histRapidityTC->fill(rapidityT, weight); 
      _histRapiditySC->fill(rapidityS, weight); 
      _histPtTInC->fill(pTinT, weight);
      _histPtTOutC->fill(pToutT, weight);
      _histPtSInC->fill(pTinS, weight);
      _histPtSOutC->fill(pToutS, weight);
    }
    const double avgPtInT  = (numParticles > 0) ? totalPtInT/fabs(numParticles) : -1;
    const double avgPtOutT = (numParticles > 0) ? totalPtOutT/fabs(numParticles) : -1;
    _histPtTInVsXp->fill(avgPtInT, weight); 
    _histPtTOutVsXp->fill(avgPtOutT, weight); 


    // Finished...
    log << Log::DEBUG << "Finished analyzing" << endl;
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

    _histPtTOutVsXp   = bookHistogram1D(9, 1, 1, "Mean out-of-plane p_T in GeV w.r.t. thrust axes vs. x_p (charged)"); // binned in Xp
    _histPtTInVsXp    = bookHistogram1D(10, 1, 1, "Mean in-plane p_T in GeV w.r.t. thrust axes vs. x_p (charged)"); // binned in Xp

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
    /// @todo Normalize most to avg number of particles (charged or otherwise)
    
    normalize(_histPtTInC);
    normalize(_histPtTInCN);
    normalize(_histPtTOutC); 
    normalize(_histPtTOutCN); 
    normalize(_histPtSInC);
    normalize(_histPtSInCN);
    normalize(_histPtSOutC); 
    normalize(_histPtSOutCN); 
    
    normalize(_histRapidityTC); 
    normalize(_histRapidityTCN); 
    normalize(_histRapiditySC); 
    normalize(_histRapiditySCN); 
    
    normalize(_histLogScaledMom); 
    normalize(_histScaledMom); 

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
