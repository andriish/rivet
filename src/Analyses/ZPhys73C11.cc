// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
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
  void ZPhys73C11::analyze(const Event & event) {
    Log& log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;

    // Analyse and print some info
    const Multiplicity& m = event.applyProjection(mult);
    log << Log::INFO << "Total charged multiplicity    = " << m.totalChargedMultiplicity() << endl;

    //Analyse the event shape info
    const Sphericity& s = event.applyProjection(spher);
    log << Log::INFO << "Sphericity    = " << s.sphericity() << endl;
    log << Log::INFO << "Aplanarity    = " << s.aplanarity() << endl;
    log << Log::INFO << "Planarity     = " << s.planarity() << endl;

    // Fill histograms here, and scale them later
    const double weight = event.weight();
    histChTot_->fill(m.totalChargedMultiplicity(), weight);
    histSphericity_->fill(s.sphericity(), weight);
    histPlanarity_->fill(s.planarity(), weight);
    histAplanarity_->fill(s.aplanarity(), weight);

    // Finished...
    log << Log::DEBUG << "Finished analyzing" << endl;
  }


  // Finalize
  void ZPhys73C11::finalize() { 
    // Normalize the histogram areas to 1
    normalize(histChTot_);
    normalize(histSphericity_);
    normalize(histPlanarity_); 
    normalize(histAplanarity_);
  }

}
