// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S7837160.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2008_S7837160::D0_2008_S7837160()
    : Analysis("D0_2008_S7837160")
  {
    // Run II W charge asymmetry
    setBeams(PROTON, ANTIPROTON);

    // Leading electrons
    FinalState fs(-5.0, 5.0);

    LeadingParticlesFinalState efs(fs);
    efs.addParticleId(ELECTRON).addParticleId(POSITRON);
    addProjection(efs, "WDecayE");

    LeadingParticlesFinalState nufs(fs);
    nufs.addParticleId(NU_E).addParticleId(NU_EBAR);
    addProjection(nufs, "WDecayNu");

    // Final state w/o electron
    VetoedFinalState vfs(fs);
    /// @todo A better way would be to have a "only photons FS". Add this projection.
    /// @todo Make it easier to exclude all neutrinos
    vfs.addVetoOnThisFinalState(efs);
    vfs.addVetoPairId(NU_E).addVetoPairId(NU_MU).addVetoPairId(NU_TAU);
    addProjection(vfs, "NoElectronFS");
  } 



  // Book histograms
  void D0_2008_S7837160::init() {
    _h_dsigplus_deta_25_35  = bookHistogram1D("dsigplus_deta_25_35", 10, 0.0, 3.2);
    _h_dsigminus_deta_25_35 = bookHistogram1D("dsigminus_deta_25_35", 10, 0.0, 3.2);
    _h_dsigplus_deta_35     = bookHistogram1D("dsigplus_deta_35", 10, 0.0, 3.2);
    _h_dsigminus_deta_35    = bookHistogram1D("dsigminus_deta_35", 10, 0.0, 3.2);
    _h_dsigplus_deta_25     = bookHistogram1D("dsigplus_deta_25", 10, 0.0, 3.2);
    _h_dsigminus_deta_25    = bookHistogram1D("dsigminus_deta_25", 10, 0.0, 3.2);
  }


  // Do the analysis 
  void D0_2008_S7837160::analyze(const Event & event) {
    const double weight = event.weight();

    // Find the W decay products
    const FinalState& efs = applyProjection<FinalState>(event, "WDecayE");
    const FinalState& nufs = applyProjection<FinalState>(event, "WDecayNu");

    // If there is no e/nu_e pair in the FinalState, skip the event
    if (efs.particles().size() < 1 || nufs.particles().size() < 1) {
      getLog() << Log::DEBUG << "No e/nu_e pair found " << endl;
      vetoEvent;
    }

    // Identify leading nu and electron
    ParticleVector es = efs.particles();
    sort(es.begin(), es.end(), cmpParticleByEt);
    Particle leading_e = es[0];
    //
    ParticleVector nus = nufs.particles();
    sort(nus.begin(), nus.end(), cmpParticleByEt);
    Particle leading_nu = nus[0];

    // Require that the neutrino has Et > 25 GeV
    const FourMomentum nu = leading_nu.momentum();
    if (nu.Et()/GeV < 25) {
      getLog() << Log::DEBUG << "Neutrino fails Et cut" << endl;
      vetoEvent;
    }

    // Get "raw" electron 4-momentum and add back in photons that could have radiated from the electron
    FourMomentum e = leading_e.momentum();
    const ParticleVector allparts = applyProjection<FinalState>(event, "NoElectronFS").particles();
    const double HALO_RADIUS = 0.2;
    foreach (const Particle& p, allparts) {
      if (p.pdgId() == PHOTON) {
        const double pho_eta = p.momentum().pseudorapidity();
        const double pho_phi = p.momentum().azimuthalAngle();
        if (deltaR(e.pseudorapidity(), e.azimuthalAngle(), pho_eta, pho_phi) < HALO_RADIUS) {
          e += p.momentum();
        }
      }
    }

    // Require that the electron has Et > 25 GeV
    if (e.Et()/GeV < 25) {
      getLog() << Log::DEBUG << "Electron fails Et cut" << endl;
      vetoEvent;
    }


    const double eta_e = fabs(e.pseudorapidity());
    const double et_e = e.Et();
    const int chg_e = PID::threeCharge(leading_e.pdgId());
    if (et_e/GeV < 35) {
      // 25 < ET < 35
      if (chg_e < 0) {
        _h_dsigminus_deta_25_35->fill(eta_e, weight);
      } else {
        _h_dsigplus_deta_25_35->fill(eta_e, weight);
      }
    } else {
      // ET > 35
      if (chg_e < 0) {
        _h_dsigminus_deta_35->fill(eta_e, weight);
      } else {
        _h_dsigplus_deta_35->fill(eta_e, weight);
      }
    }
    // Inclusive: ET > 25
    if (chg_e < 0) {
      _h_dsigminus_deta_25->fill(eta_e, weight);
    } else {
      _h_dsigplus_deta_25->fill(eta_e, weight);
    }
  }


  // Finalize
  void D0_2008_S7837160::finalize() {
    // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta) for each Et region
    AIDA::IHistogramFactory& hf = histogramFactory();

    const string basetitle = "W charge asymmetry for ";
    const string xlabel = "$|\\eta|$ of leading electron";
    const string ylabel = "A = "
      "$(\\frac{\\mathrm{d}{\\sigma^+}}{\\mathrm{d}{|\\eta|}} - \\frac{\\mathrm{d}{\\sigma^-}}{\\mathrm{d}{|\\eta|}}) / "
      "(\\frac{\\mathrm{d}{\\sigma^+}}{\\mathrm{d}{|\\eta|}} + \\frac{\\mathrm{d}{\\sigma^-}}{\\mathrm{d}{|\\eta|}})$";

    IHistogram1D* num25_35 = hf.subtract("/num25_35", *_h_dsigplus_deta_25_35, *_h_dsigminus_deta_25_35);
    IHistogram1D* denom25_35 = hf.add("/denom25_35", *_h_dsigplus_deta_25_35, *_h_dsigminus_deta_25_35);
    assert(num25_35 && denom25_35);
    IDataPointSet* tot25_35 = hf.divide(histoDir() + "/d01-x01-y01", *num25_35, *denom25_35);
    tot25_35->setTitle(basetitle + "$25 < E_\\perp < 35$ GeV");
    tot25_35->setXTitle(xlabel);
    tot25_35->setYTitle(ylabel);
    hf.destroy(num25_35);
    hf.destroy(denom25_35);
    //
    IHistogram1D* num35 = hf.subtract("/num35", *_h_dsigplus_deta_35, *_h_dsigminus_deta_35);
    IHistogram1D* denom35 = hf.add("/denom35", *_h_dsigplus_deta_35, *_h_dsigminus_deta_35);
    assert(num35 && denom35);
    IDataPointSet* tot35 = hf.divide(histoDir() + "/d02-x01-y01", *num35, *denom35);
    tot35->setTitle(basetitle + "$E_\\perp > 35$ GeV");
    tot35->setXTitle(xlabel);
    tot35->setYTitle(ylabel);
    hf.destroy(num35);
    hf.destroy(denom35);
    //
    IHistogram1D* num25 = hf.subtract("/num25", *_h_dsigplus_deta_25, *_h_dsigminus_deta_25);
    IHistogram1D* denom25 = hf.add("/denom25", *_h_dsigplus_deta_25, *_h_dsigminus_deta_25);
    assert(num25 && denom25);
    IDataPointSet* tot25 = hf.divide(histoDir() + "/d03-x01-y01", *num25, *denom25);
    tot25->setTitle(basetitle + "$E_\\perp > 35$ GeV");
    tot25->setXTitle(xlabel);
    tot25->setYTitle(ylabel);
    hf.destroy(num25);
    hf.destroy(denom25);
  }


}
