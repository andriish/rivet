#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Analyses/MC_LHC_DIJET.hh"

namespace Rivet {

  /// Default constructor
  MC_LHC_DIJET::MC_LHC_DIJET() 
  : Analysis("MC_LHC_DIJET") {
    FinalState fs;
    ChargedFinalState cfs;
    addProjection(fs, "FS");
    addProjection(cfs, "CFS");
    addProjection(FastJets(fs, FastJets::KT, 0.7), "Jets");
    addProjection(FastJets(cfs, FastJets::KT, 0.7), "ChargedJets");
  }


  /// @name Analysis methods
  //@{
  void MC_LHC_DIJET::init() { 
    _hist_jetcount = bookHistogram1D("d01-x01-y01", 5, 0., 10.);
    _hist_jetpt = bookHistogram1D("d02-x01-y01", 30, 30.,100.);
    _hist_jetptlog = bookHistogram1D("d03-x01-y01", 20, 0.,8.);
    _hist_leadingjetpt = bookHistogram1D("d04-x01-y01", 25, 30.,100.);
    _hist_secondleadingjetpt = bookHistogram1D("d05-x01-y01", 25, 30.,100.);
    _hist_jetphi = bookHistogram1D("d06-x01-y01",24, 0., 6.4);
    _hist_jeteta = bookHistogram1D("d07-x01-y01", 30, -6., 6.);
    _hist_jetdphi = bookHistogram1D("d08-x01-y01", 24, 0., 6.4);
    _hist_jetdeta = bookHistogram1D("d09-x01-y01", 24, 0., 6.);
    _hist_chargemultiplicity = bookHistogram1D("d10-x01-y01",30, 0.5, 250.5);
    _hist_chargemeanpt = bookHistogram1D("d11-x01-y01", 25, 0., 10.);    
    _hist_chargept = bookHistogram1D("d12-x01-y01", 32, 0., 25.);
    _hist_chargelogpt = bookHistogram1D("d13-x01-y01", 32, 0., 6.);    
    _hist_chargermspt = bookHistogram1D("d14-x01-y01", 32, 0., 10.);
  }



  void MC_LHC_DIJET::analyze(const Event& event) {
    const FastJets& fastjets = applyProjection<FastJets>(event, "Jets");
    const Jets jets = fastjets.jetsByPt(20.);
    const double weight = event.weight();
    
    if (jets.size() < 2 || jets.size() >= 3) vetoEvent;
    const double angle = fabs(jets[1].momentum().azimuthalAngle() - jets[0].momentum().azimuthalAngle());
    const double prapidity = fabs(jets[1].momentum().pseudorapidity() - jets[0].momentum().pseudorapidity());
    _hist_jetcount->fill(jets.size(), weight);
    _hist_leadingjetpt->fill(jets[0].momentum().pT(), weight);
    _hist_secondleadingjetpt->fill(jets[1].momentum().pT(), weight);
    _hist_jetdphi->fill(angle , weight);
    _hist_jetdeta->fill(prapidity, weight);
    
    foreach(Jet j, fastjets.jetsByPt(20*GeV)) {
      _hist_jetpt->fill(j.momentum().pT(), weight);
      _hist_jetptlog->fill(log(j.momentum().pT()), weight);
      _hist_jetphi->fill(j.momentum().azimuthalAngle(), weight);
      _hist_jeteta->fill(j.momentum().pseudorapidity(), weight);	
    }
    
    const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
    // const FastJets& cfastjets = applyProjection<FastJets>(event, "ChargedJets");
    // const Jets cjets = cfastjets.jetsByPt(20.);
    double meanpt = 0.;
    double rmspt = 0.;
    // foreach(Jet cj, cfastjets.jetsByPt(20.)){
    _hist_chargemultiplicity->fill(cfs.particles().size(), weight);
    foreach(Particle cp, cfs.particles()) {
      meanpt= meanpt + cp.momentum().pT();
      rmspt = rmspt + (cp.momentum().pT()*cp.momentum().pT());
      _hist_chargept->fill(cp.momentum().pT(), weight);
      _hist_chargelogpt->fill(log(cp.momentum().pT()), weight);
    }
    meanpt = meanpt / cfs.particles().size();
    _hist_chargemeanpt->fill(meanpt, weight);
    rmspt = sqrt(rmspt / cfs.particles().size());
    _hist_chargermspt->fill(rmspt, weight);
    // }
  }
  

  void MC_LHC_DIJET::finalize() {
    /// @todo Normalise!
  }

  //@}



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_LHC_DIJET> plugin_MC_LHC_DIJET;

}
