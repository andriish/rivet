// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/RivetCLHEP.hh"
#include "Rivet/Analyses/ZPhys73C11.hh"
#include "HepPDT/ParticleID.hh"

using namespace AIDA;
using namespace HepMC;


namespace Rivet {

  void ZPhys73C11::init() {
    // Inclusive single-particle distributions
    _histLogScaledMom = bookHistogram1D("IncSingleLogScaledMom", "Inclusive single particle log(1/x_p)", 60, 0.0, 6.0);
    _histScaledMom    = bookHistogram1D("IncSingleScaledMom", "Inclusive single particle x_p = |p|/|p_beam|", 50, 0.0, 1.0);
    _histRapidityT    = bookHistogram1D("IncSingleRapidityT", "Inclusive single particle rapidity w.r.t. thrust axis, y_T", 30, 0.0, 6.0);
    _histRapidityS    = bookHistogram1D("IncSingleRapidityS", "Inclusive single particle rapidity w.r.t. sphericity axis, y_S", 30, 0.0, 6.0);
    _histPtTIn        = bookHistogram1D("IncSinglePtTIn", "Inclusive single particle in-plane p_T w.r.t. thrust axis (GeV)", 50, 0.0, 14.0);
    _histPtTOut       = bookHistogram1D("IncSinglePtTOut", "Inclusive single particle out-of-plane p_T w.r.t. thrust-minor axis (GeV)", 50, 0.0, 3.5);
    _histPtTInVsXp    = bookHistogram1D("IncSinglePtTInMeanVsXp", "Inclusive single particle in-plane p_T (GeV) vs. x_p", 50, 0.0, 3.5);
    _histPtTOutVsXp   = bookHistogram1D("IncSinglePtTOutMeanVsXp", "Inclusive single particle out-of-plane p_T (GeV) vs. x_p", 50, 0.0, 3.5);

    // Event shape distributions
    _hist1MinusT      = bookHistogram1D("ESOneMinusThrust", "1-thrust, 1-T", 50, 0.0, 0.5);
    _histTMajor       = bookHistogram1D("ESMajor", "Thrust major, M", 50, 0.0, 0.65);
    _histTMinor       = bookHistogram1D("ESMinor", "Thrust minor, m", 50, 0.0, 0.4);
    _histOblateness   = bookHistogram1D("ESOblateness", "Oblateness = M - m", 50, 0.0, 0.5);

    _histDiffRate2Durham = bookHistogram1D("DiffRate2Durham", "Differential 2-jet rate with Durham algorithm, D_2^Durham", 20, 0.0, 0.3); // binned in y_cut
    _histDiffRate2Jade   = bookHistogram1D("DiffRate2Jade", "Differential 2-jet rate with Jade algorithm, D_2^Jade", 20, 0.0, 0.3); // binned in y_cut
    _histDiffRate2JadeCN = bookHistogram1D("DiffRate2JadeCN", "Differential 2-jet rate with Jade algorithm from charged and neutral particles, D_2,CN^Jade", 20, 0.0, 0.3); // binned in y_cut
    _histDiffRate3Durham = bookHistogram1D("DiffRate3Durham", "Differential 3-jet rate with Durham algorithm, D_3^Durham", 20, 0.0, 0.1); // binned in y_cut
    _histDiffRate3Jade   = bookHistogram1D("DiffRate3Jade", "Differential 3-jet rate with Jade algorithm from charged and neutral particles, D_3,CN^Jade", 20, 0.0, 0.1); // binned in y_cut
    _histDiffRate3JadeCN = bookHistogram1D("DiffRate3JadeCN", "Differential 3-jet rate with Jade algorithm, D_3^Jade", 20, 0.0, 0.1); // binned in y_cut
    _histDiffRate4Durham = bookHistogram1D("DiffRate4Durham", "Differential 4-jet rate with Durham algorithm, D_4^Durham", 20, 0.0, 0.03); // binned in y_cut
    _histDiffRate4Jade   = bookHistogram1D("DiffRate4Jade", "Differential 4-jet rate with Jade algorithm, D_4^Jade", 20, 0.0, 0.05); // binned in y_cut

    _histSphericity  = bookHistogram1D("ESSphericity", "Sphericity, S", 20, 0.0, 0.8);
    _histAplanarity  = bookHistogram1D("ESAplanarity", "Aplanarity, A", 20, 0.0, 0.3);
    _histPlanarity   = bookHistogram1D("ESPlanarity", "Planarity, P", 20, 0.0, 0.5);
    _histSphericityCN  = bookHistogram1D("ESSphericityCN", "Sphericity from charged and neutral particles, S_CN", 20, 0.0, 0.8);
    _histAplanarityCN  = bookHistogram1D("ESAplanarityCN", "Aplanarity from charged and neutral particles, A_CN", 20, 0.0, 0.3);

    _histHemiMassD = bookHistogram1D("HemiMassDiff", "Difference in hemisphere masses, M_d^2/E_vis^2", 20, 0.0, 0.4);
    _histHemiMassH = bookHistogram1D("HemiMassHeavy", "Heavy hemisphere masses, M_h^2/E_vis^2", 20, 0.0, 0.4);
    _histHemiMassL = bookHistogram1D("HemiMassLight", "Light hemisphere masses, M_l^2/E_vis^2", 20, 0.0, 0.12);

    _histHemiBroadW = bookHistogram1D("HemiBroadW", "Wide hemisphere broadening, B_max", 20, 0.0, 0.3);
    _histHemiBroadN = bookHistogram1D("HemiBroadN", "Narrow hemisphere broadening, B_min", 20, 0.0, 0.18);
    _histHemiBroadT = bookHistogram1D("HemiBroadT", "Total hemisphere broadening, B_sum", 20, 0.0, 0.35);
    _histHemiBroadD = bookHistogram1D("HemiBroadD", "Difference in hemisphere broadening, B_diff", 20, 0.0, 0.30);

    _histCParam = bookHistogram1D("CParam", "C parameter, C", 20, 0.0, 0.9);
    _histDParam = bookHistogram1D("DParam", "D parameter, D", 20, 0.0, 0.8);

    _histEEC = bookHistogram1D("EEC", "Energy-energy correlation, EEC", 40, -1.0, 1.0); // binned in cos(chi)
    _histAEEC = bookHistogram1D("AEEC", "Asymmetry of the energy-energy correlation, AEEC", 20, -1.0, 0.0); // binned in cos(chi)


    // Identified particle distributions
    // ?

  }


  // Do the analysis
  void ZPhys73C11::analyze(const Event& e) {
    Log& log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;

    // Get beams and average beam momentum
    /// @todo Some problem with how the beams projection's destructor handles particles: segfaults.
    const ParticlePair& beams = e.applyProjection(_beamsproj)();
    const double meanBeamMom = ( beams.first.getMomentum().vect().mag() + beams.second.getMomentum().vect().mag() ) / 2.0;

    // Calculate event-wise shape distributions
    const Thrust& thrust = e.applyProjection(_thrustproj);
    const Sphericity& sphericityC = e.applyProjection(_cspherproj);
    const Sphericity& sphericityCN = e.applyProjection(_cnspherproj);
    const ParisiTensor& parisi = e.applyProjection(_parisiproj);

    // Get event weight for histo filling
    const double weight = e.weight();

    // Iterate over all the final state particles.
    const FinalState& cfs = e.applyProjection(_cfsproj);
    for (ParticleVector::const_iterator p = cfs.particles().begin(); p != cfs.particles().end(); ++p) {
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
      _histPtTInVsXp->fill(pTinS, weight); 
      _histPtTOutVsXp->fill(pToutS, weight); 
    }

    // Fill histograms.
    _hist1MinusT->fill(1 - thrust.thrust(), weight); 
    _histTMajor->fill(thrust.thrustMajor(), weight); 
    _histTMinor->fill(thrust.thrustMinor(), weight); 
    _histOblateness->fill(thrust.oblateness(), weight);

    //_histDiffRate2Durham->fill(, weight); 
    //_histDiffRate2Jade->fill(, weight); 
    //_histDiffRate2JadeCN->fill(, weight); 
    //_histDiffRate3Durham->fill(, weight);
    //_histDiffRate3Jade->fill(, weight); 
    //_histDiffRate3JadeCN->fill(, weight);
    //_histDiffRate4Durham->fill(, weight);
    //_histDiffRate4Jade->fill(, weight); 
    
    _histSphericity->fill(sphericityC.sphericity(), weight); 
    _histAplanarity->fill(sphericityC.aplanarity(), weight); 
    _histPlanarity->fill(sphericityC.planarity(), weight); 
    _histSphericityCN->fill(sphericityCN.sphericity(), weight); 
    _histAplanarityCN->fill(sphericityCN.aplanarity(), weight); 

    //_histHemiMassD->fill(, weight); 
    //_histHemiMassH->fill(, weight); 
    //_histHemiMassL->fill(, weight); 

    //_histHemiBroadW->fill(, weight); 
    //_histHemiBroadN->fill(, weight); 
    //_histHemiBroadT->fill(, weight); 
    //_histHemiBroadD->fill(, weight); 
    
    _histCParam->fill(parisi.C(), weight); 
    _histDParam->fill(parisi.D(), weight); 
    
    //_histEEC->fill(, weight); 
    //_histAEEC->fill(, weight); 

    // Finished...
    log << Log::DEBUG << "Finished analyzing" << endl;
  }


  // Finalize
  void ZPhys73C11::finalize() { 
    // Normalize the histogram areas to 1
//     normalize(_histLogScaledMom); 
//     normalize(_histScaledMom); 
//     normalize(_histRapidityT); 
//     normalize(_histRapidityS); 
//     normalize(_histPtTIn); 
//     normalize(_histPtTOut); 
//     normalize(_histPtTInVsXp); 
//     normalize(_histPtTOutVsXp); 
    
//     normalize(_hist1MinusT); 
//     normalize(_histTMajor); 
//     normalize(_histTMinor); 
//     normalize(_histOblateness); 
    
//     normalize(_histDiffRate2Durham); 
//     normalize(_histDiffRate2Jade); 
//     normalize(_histDiffRate2JadeCN);
//     normalize(_histDiffRate3Durham);
//     normalize(_histDiffRate3Jade); 
//     normalize(_histDiffRate3JadeCN);
//     normalize(_histDiffRate4Durham);
//     normalize(_histDiffRate4Jade); 
    
//     normalize(_histSphericity); 
//     normalize(_histAplanarity); 
//     normalize(_histPlanarity); 
//     normalize(_histSphericityCN); 
//     normalize(_histAplanarityCN); 
    
//     normalize(_histHemiMassD); 
//     normalize(_histHemiMassH); 
//     normalize(_histHemiMassL); 
    
//     normalize(_histHemiBroadW); 
//     normalize(_histHemiBroadN); 
//     normalize(_histHemiBroadT); 
//     normalize(_histHemiBroadD); 
    
//     normalize(_histCParam); 
//     normalize(_histDParam); 
    
//     normalize(_histEEC); 
//     normalize(_histAEEC); 
  }

}
