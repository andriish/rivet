// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/RivetCLHEP.hh"
#include "Rivet/Analyses/ZPhys73C11.hh"
#include "HepPDT/ParticleID.hh"

using namespace AIDA;
using namespace HepMC;


namespace Rivet {


  void ZPhys73C11::analyze(const Event& e) {
    Log& log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;

    // Get beams and average beam momentum
    const ParticlePair& beams = e.applyProjection(_beamsproj).getBeams();
    const double meanBeamMom = ( beams.first.getMomentum().vect().mag() + 
                                 beams.second.getMomentum().vect().mag() ) / 2.0;

    // Get event weight for histo filling
    const double weight = e.weight();

    const Thrust& thrust = e.applyProjection(_thrustproj);
    _hist1MinusT->fill(1 - thrust.thrust(), weight); 
    _histTMajor->fill(thrust.thrustMajor(), weight); 
    _histTMinor->fill(thrust.thrustMinor(), weight); 
    _histOblateness->fill(thrust.oblateness(), weight);

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

    const Sphericity& sphericityC = e.applyProjection(_cspherproj);
    _histSphericity->fill(sphericityC.sphericity(), weight); 
    _histAplanarity->fill(sphericityC.aplanarity(), weight); 
    _histPlanarity->fill(sphericityC.planarity(), weight); 

    const Sphericity& sphericityCN = e.applyProjection(_cnspherproj);
    _histSphericityCN->fill(sphericityCN.sphericity(), weight); 
    _histAplanarityCN->fill(sphericityCN.aplanarity(), weight); 

    const ParisiTensor& parisi = e.applyProjection(_parisiproj);    
    _histCParam->fill(parisi.C(), weight); 
    _histDParam->fill(parisi.D(), weight); 

    const Hemispheres& hemi = e.applyProjection(_hemiproj);
    _histHemiMassH->fill(hemi.getScaledM2high(), weight); 
    _histHemiMassL->fill(hemi.getScaledM2low(), weight); 
    _histHemiMassD->fill(hemi.getScaledM2diff(), weight); 
    _histHemiBroadW->fill(hemi.getBmax(), weight); 
    _histHemiBroadN->fill(hemi.getBmin(), weight); 
    _histHemiBroadT->fill(hemi.getBsum(), weight); 
    _histHemiBroadD->fill(hemi.getBdiff(), weight); 

    //const EEC& eec = e.applyProjection(_eecproj);
    //_histEEC->fill(, weight); 
    //_histAEEC->fill(, weight); 


    // Iterate over all the final state particles.
    double totalPtInT(0), totalPtOutT(0);
    const FinalState& cfs = e.applyProjection(_cfsproj);
    const size_t numParticles = cfs.particles().size();
    for (ParticleVector::const_iterator p = cfs.particles().begin(); p != cfs.particles().end(); ++p) {
      /// @todo Add missing CN variants

      // Get momentum and energy of each particle.
      const Vector3 mom3 = p->getMomentum().vect();
      const double mom = mom3.mag();
      const double energy = p->getMomentum().getT();

      // Calculate scaled momenta.
      const double scaledMom = mom/meanBeamMom;
      const double logInvScaledMom = -log10(scaledMom);

      // Get momenta components w.r.t. thrust and sphericity.
      const double momT = thrust.thrustAxis().dot(mom3);
      const double momS = sphericityC.sphericityAxis().dot(mom3);
      const double pTinT = mom3.dot(thrust.thrustMajorAxis());
      const double pToutT = mom3.dot(thrust.thrustMinorAxis());
      const double pTinS = mom3.dot(sphericityC.sphericityMinorAxis());
      const double pToutS = mom3.dot(sphericityC.sphericityMinorAxis());
      totalPtInT += pTinT; 
      totalPtOutT += pToutT;

      // Calculate rapidities w.r.t. thrust and sphericity.
      const double rapidityT = 0.5 * (energy + momT) / (energy - momT);
      const double rapidityS = 0.5 * (energy + momS) / (energy - momS);

      // Fill histograms.
      _histLogScaledMom->fill(logInvScaledMom, weight); 
      _histScaledMom->fill(scaledMom, weight); 
      _histRapidityT->fill(rapidityT, weight); 
      _histRapidityS->fill(rapidityS, weight); 
      _histPtTIn->fill(pTinT, weight);
      _histPtTOut->fill(pToutT, weight);
      _histPtSIn->fill(pTinS, weight);
      _histPtSOut->fill(pToutS, weight);
    }
    const double avgPtInT  = (numParticles > 0) ? totalPtInT/fabs(numParticles) : -1;
    const double avgPtOutT = (numParticles > 0) ? totalPtOutT/fabs(numParticles) : -1;
    _histPtTInVsXp->fill(avgPtInT, weight); 
    _histPtTOutVsXp->fill(avgPtOutT, weight); 


    // Finished...
    log << Log::DEBUG << "Finished analyzing" << endl;
  }



  void ZPhys73C11::init() {
    _histPtTIn        = bookHistogram1D(1, 1, 1, "In-plane p_T in GeV w.r.t. thrust axes (charged)");
    _histPtTInCN      = bookHistogram1D(1, 1, 2, "In-plane p_T in GeV w.r.t. thrust axes (charged and neutral)");
    _histPtTOut       = bookHistogram1D(2, 1, 1, "Out-of-plane p_T in GeV w.r.t. thrust axes (charged)");
    _histPtTOutCN     = bookHistogram1D(2, 1, 2, "Out-of-plane p_T in GeV w.r.t. thrust axes (charged and neutral)");
    _histPtSIn        = bookHistogram1D(3, 1, 1, "In-plane p_T in GeV w.r.t. sphericity axes (charged)");
    _histPtSInCN      = bookHistogram1D(3, 1, 2, "In-plane p_T in GeV w.r.t. sphericity axes (charged and neutral)");
    _histPtSOut       = bookHistogram1D(4, 1, 1, "Out-of-plane p_T in GeV w.r.t. sphericity axes (charged)");
    _histPtSOutCN     = bookHistogram1D(4, 1, 2, "Out-of-plane p_T in GeV w.r.t. sphericity axes (charged and neutral)");

    _histRapidityT    = bookHistogram1D(5, 1, 1, "Rapidity w.r.t. thrust axes, y_T (charged)");
    _histRapidityTCN  = bookHistogram1D(5, 1, 2, "Rapidity w.r.t. thrust axes, y_T (charged and neutral)");
    _histRapidityS    = bookHistogram1D(6, 1, 1, "Rapidity w.r.t. sphericity axes, y_S (charged)");
    _histRapiditySCN  = bookHistogram1D(6, 1, 2, "Rapidity w.r.t. sphericity axes, y_S (charged and neutral)");

    _histScaledMom    = bookHistogram1D(7, 1, 1, "Scaled momentum, x_p = |p|/|p_beam| (charged)");
    _histLogScaledMom = bookHistogram1D(8, 1, 1, "Log of scaled momentum, log(1/x_p) (charged)");

    _histPtTOutVsXp   = bookHistogram1D(9, 1, 1, "Mean out-of-plane p_T in GeV w.r.t. thrust axes vs. x_p (charged)"); // binned in Xp
    _histPtTInVsXp    = bookHistogram1D(10, 1, 1, "Mean in-plane p_T in GeV w.r.t. thrust axes vs. x_p (charged)"); // binned in Xp

    _hist1MinusT      = bookHistogram1D(11, 1, 1, "1-thrust, 1-T (charged)");
    _hist1MinusTCN    = bookHistogram1D(11, 1, 2, "1-thrust, 1-T (charged and neutral)");
    _histTMajor       = bookHistogram1D(12, 1, 1, "Thrust major, M (charged)");
    _histTMajorCN     = bookHistogram1D(12, 1, 2, "Thrust major, M (charged and neutral)");
    _histTMinor       = bookHistogram1D(13, 1, 1, "Thrust minor, m (charged)");
    _histTMinorCN     = bookHistogram1D(13, 1, 2, "Thrust minor, m (charged and neutral)");
    _histOblateness   = bookHistogram1D(14, 1, 1, "Oblateness = M - m (charged)");
    _histOblatenessCN = bookHistogram1D(14, 1, 2, "Oblateness = M - m (charged and neutral)");

    _histSphericity   = bookHistogram1D(15, 1, 1, "Sphericity, S (charged)");
    _histSphericityCN = bookHistogram1D(15, 1, 2, "Sphericity, S (charged and neutral)");
    _histAplanarity   = bookHistogram1D(16, 1, 1, "Aplanarity, A (charged)");
    _histAplanarityCN = bookHistogram1D(16, 1, 2, "Aplanarity, A (charged and neutral)");
    _histPlanarity    = bookHistogram1D(17, 1, 1, "Planarity, P (charged)");
    _histPlanarityCN  = bookHistogram1D(17, 1, 2, "Planarity, P (charged and neutral)");

    _histCParam       = bookHistogram1D(18, 1, 1, "C parameter (charged)");
    _histCParamCN     = bookHistogram1D(18, 1, 2, "C parameter (charged and neutral)");
    _histDParam       = bookHistogram1D(19, 1, 1, "D parameter (charged)");
    _histDParamCN     = bookHistogram1D(19, 1, 2, "D parameter (charged and neutral)");

    _histHemiMassH    = bookHistogram1D(20, 1, 1, "Heavy hemisphere masses, M_h^2/E_vis^2 (charged)");
    _histHemiMassHCN  = bookHistogram1D(20, 1, 2, "Heavy hemisphere masses, M_h^2/E_vis^2 (charged and neutral)");
    _histHemiMassL    = bookHistogram1D(21, 1, 1, "Light hemisphere masses, M_l^2/E_vis^2 (charged)");
    _histHemiMassLCN  = bookHistogram1D(21, 1, 2, "Light hemisphere masses, M_l^2/E_vis^2 (charged and neutral)");
    _histHemiMassD    = bookHistogram1D(22, 1, 1, "Difference in hemisphere masses, M_d^2/E_vis^2 (charged)");
    _histHemiMassDCN  = bookHistogram1D(22, 1, 2, "Difference in hemisphere masses, M_d^2/E_vis^2 (charged and neutral)");

    _histHemiBroadW   = bookHistogram1D(23, 1, 1, "Wide hemisphere broadening, B_max (charged)");
    _histHemiBroadWCN = bookHistogram1D(23, 1, 2, "Wide hemisphere broadening, B_max (charged and neutral)");
    _histHemiBroadN   = bookHistogram1D(24, 1, 1, "Narrow hemisphere broadening, B_min (charged)");
    _histHemiBroadNCN = bookHistogram1D(24, 1, 2, "Narrow hemisphere broadening, B_min (charged and neutral)");
    _histHemiBroadT   = bookHistogram1D(25, 1, 1, "Total hemisphere broadening, B_sum (charged)");
    _histHemiBroadTCN = bookHistogram1D(25, 1, 2, "Total hemisphere broadening, B_sum (charged and neutral)");
    _histHemiBroadD   = bookHistogram1D(26, 1, 1, "Difference in hemisphere broadening, B_diff (charged)");
    _histHemiBroadDCN = bookHistogram1D(26, 1, 2, "Difference in hemisphere broadening, B_diff (charged and neutral)");

    _histDiffRate2Durham   = bookHistogram1D(27, 1, 1, "Differential 2-jet rate with Durham algorithm, D_2^Durham (charged)"); // binned in y_cut
    _histDiffRate2DurhamCN = bookHistogram1D(27, 1, 2, "Differential 2-jet rate with Durham algorithm, D_2^Durham (charged and neutral)"); // binned in y_cut
    _histDiffRate2Jade     = bookHistogram1D(28, 1, 1, "Differential 2-jet rate with Jade algorithm, D_2^Jade (charged)"); // binned in y_cut
    _histDiffRate2JadeCN   = bookHistogram1D(28, 1, 2, "Differential 2-jet rate with Jade algorithm, D_2^Jade (charged and neutral)"); // binned in y_cut
    _histDiffRate3Durham   = bookHistogram1D(29, 1, 1, "Differential 3-jet rate with Durham algorithm, D_3^Durham (charged)"); // binned in y_cut
    _histDiffRate3DurhamCN = bookHistogram1D(29, 1, 2, "Differential 3-jet rate with Durham algorithm, D_3^Durham (charged and neutral)"); // binned in y_cut
    _histDiffRate3Jade     = bookHistogram1D(30, 1, 1, "Differential 3-jet rate with Jade algorithm, D_3^Jade (charged)"); // binned in y_cut
    _histDiffRate3JadeCN   = bookHistogram1D(30, 1, 2, "Differential 3-jet rate with Jade algorithm, D_3^Jade (charged and neutral)"); // binned in y_cut
    _histDiffRate4Durham   = bookHistogram1D(31, 1, 1, "Differential 4-jet rate with Durham algorithm, D_4^Durham (charged)"); // binned in y_cut
    _histDiffRate4DurhamCN = bookHistogram1D(31, 1, 2, "Differential 4-jet rate with Durham algorithm, D_4^Durham (charged and neutral)"); // binned in y_cut
    _histDiffRate4Jade     = bookHistogram1D(32, 1, 1, "Differential 4-jet rate with Jade algorithm, D_4^Jade (charged)"); // binned in y_cut
    _histDiffRate4JadeCN   = bookHistogram1D(32, 1, 2, "Differential 4-jet rate with Jade algorithm, D_4^Jade (charged and neutral)"); // binned in y_cut

    _histEEC  = bookHistogram1D(33, 1, 1, "Energy-energy correlation, EEC (charged)"); // binned in cos(chi)
    _histAEEC = bookHistogram1D(34, 1, 1, "Asymmetry of the energy-energy correlation, AEEC (charged)"); // binned in cos(chi)


    // Identified particle distributions
    // ?
  }



  // Finalize
  void ZPhys73C11::finalize() { 
     normalize(_histPtTIn);
     normalize(_histPtTInCN);
     normalize(_histPtTOut); 
     normalize(_histPtTOutCN); 
     normalize(_histPtSIn);
     normalize(_histPtSInCN);
     normalize(_histPtSOut); 
     normalize(_histPtSOutCN); 

     normalize(_histRapidityT); 
     normalize(_histRapidityTCN); 
     normalize(_histRapidityS); 
     normalize(_histRapiditySCN); 

     normalize(_histLogScaledMom); 
     normalize(_histScaledMom); 

     normalize(_hist1MinusT); 
     normalize(_hist1MinusTCN); 
     normalize(_histTMajor); 
     normalize(_histTMajorCN); 
     normalize(_histTMinor); 
     normalize(_histTMinorCN); 
     normalize(_histOblateness); 
     normalize(_histOblatenessCN); 

    normalize(_histSphericity); 
    normalize(_histSphericityCN); 
    normalize(_histAplanarity); 
    normalize(_histAplanarityCN); 
    normalize(_histPlanarity); 
    normalize(_histPlanarityCN); 

    normalize(_histHemiMassD); 
    normalize(_histHemiMassDCN); 
    normalize(_histHemiMassH); 
    normalize(_histHemiMassHCN); 
    normalize(_histHemiMassL); 
    normalize(_histHemiMassLCN); 
    
    normalize(_histHemiBroadW); 
    normalize(_histHemiBroadWCN); 
    normalize(_histHemiBroadN); 
    normalize(_histHemiBroadNCN); 
    normalize(_histHemiBroadT); 
    normalize(_histHemiBroadTCN); 
    normalize(_histHemiBroadD); 
    normalize(_histHemiBroadDCN); 
    
    normalize(_histCParam); 
    normalize(_histCParamCN); 
    normalize(_histDParam); 
    normalize(_histDParamCN); 
    
    normalize(_histDiffRate2Durham); 
    normalize(_histDiffRate2DurhamCN); 
    normalize(_histDiffRate2Jade); 
    normalize(_histDiffRate2JadeCN);
    normalize(_histDiffRate3Durham);
    normalize(_histDiffRate3DurhamCN);
    normalize(_histDiffRate3Jade); 
    normalize(_histDiffRate3JadeCN);
    normalize(_histDiffRate4Durham);
    normalize(_histDiffRate4DurhamCN);
    normalize(_histDiffRate4Jade); 
    normalize(_histDiffRate4JadeCN); 
  }


}
