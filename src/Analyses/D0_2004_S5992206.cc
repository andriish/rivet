// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/D0_2004_S5992206.hh"


namespace Rivet {


  // Book histograms
  void D0_2004_S5992206::init() {
    _histJetAzimuth_pTmax75_100  = 
      bookHistogram1D(1, 2, 1, "Jet--jet azimuthal angle, $p_\\perp^\\text{max} \\in [75,100]\\text{ GeV}$");
    _histJetAzimuth_pTmax100_130 = 
      bookHistogram1D(2, 2, 1, "Jet--jet azimuthal angle, $p_\\perp^\\text{max} \\in [100..130]\\text{ GeV}$");
    _histJetAzimuth_pTmax130_180 = 
      bookHistogram1D(3, 2, 1, "Jet--jet azimuthal angle, $p_\\perp^\\text{max} \\in [130..180]\\text{ GeV}$");
    _histJetAzimuth_pTmax180_    = 
      bookHistogram1D(4, 2, 1, "Jet--jet azimuthal angle, $p_\\perp^\\text{max} > 180\\text{ GeV}$");
  }


  // Do the analysis
  void D0_2004_S5992206::analyze(const Event & event) {
    Log& log = getLog();

    // Analyse and print some info  
    const D0ILConeJets& jetpro = applyProjection<D0ILConeJets>(event, "Jets");
    log << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.getNJets() << endl;

    // Find vertex and check  that its z-component is < 50 cm from the nominal IP
    //const PVertex& pv = event.applyProjection(_vertexproj);
    /// @todo SEGV: either the HepMC event record is not filled properly or the F77-Wrapper functions are faulty
    // if (fabs(pv.getPrimaryVertex().position().z())/mm > 500.0) {
    //   vetoEvent(event);
    // }

    const list<FourMomentum>& jets = jetpro.getLorentzJets();
    list<FourMomentum>::const_iterator jetpTmax = jets.end();
    list<FourMomentum>::const_iterator jet2ndpTmax = jets.end();
    log << Log::DEBUG << "jetlist size = " << jets.size() << endl;
    
    size_t Njet = 0;
    for (list<FourMomentum>::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      log << Log::DEBUG << "List item pT = " << jt->pT() << " E=" << jt->E() 
          << " pz=" << jt->pz() << endl;
      if (jt->pT() > 40.0) ++Njet;
      log << Log::DEBUG << "Jet pT =" << jt->pT() << " y=" << jt->rapidity() 
          << " phi=" << jt->azimuthalAngle() << endl; 
      if (jetpTmax == jets.end() || jt->pT() > jetpTmax->pT()) {
        jet2ndpTmax = jetpTmax;
        jetpTmax = jt;
      } else if (jet2ndpTmax == jets.end() || jt->pT() > jet2ndpTmax->pT()) {
        jet2ndpTmax = jt;
      }
    }
    if (Njet >= 2) {
      log << Log::DEBUG << "Jet multiplicity after pT > 40GeV cut = " << Njet << endl; 
    }

    /// @todo Use cut constants and register these cuts
    if (jets.size() < 2 || jet2ndpTmax->pT()/GeV < 40.0) {
      vetoEvent(event);
    }
    /// @todo Use cut constants and register these cuts
    if (fabs(jetpTmax->rapidity()) > 0.5 || fabs(jet2ndpTmax->rapidity()) > 0.5) {
      vetoEvent(event);
    }
    log << Log::DEBUG << "Jet eta and pT requirements fulfilled" << endl;

    const TotalVisibleMomentum& caloMissEt = applyProjection<TotalVisibleMomentum>(event, "CalMET");
    log << Log::DEBUG << "Missing Et = " << caloMissEt.getMomentum().pT() << endl;
    /// @todo Use cut constants and register these cuts
    if (caloMissEt.getMomentum().pT() > 0.7*jetpTmax->pT()) {
      vetoEvent(event);
    }
    
    const double dphi = delta_phi(jetpTmax->azimuthalAngle(), jet2ndpTmax->azimuthalAngle());
    if (jetpTmax->pT()/GeV > 75.0 && jetpTmax->pT()/GeV <= 100.0)
      _histJetAzimuth_pTmax75_100->fill(dphi, event.weight());
    else if (jetpTmax->pT()/GeV > 100.0 && jetpTmax->pT()/GeV <= 130.0)
      _histJetAzimuth_pTmax100_130->fill(dphi, event.weight());
    else if (jetpTmax->pT()/GeV > 130.0 && jetpTmax->pT()/GeV <= 180.0)
      _histJetAzimuth_pTmax130_180->fill(dphi, event.weight());
    else if (jetpTmax->pT()/GeV > 180.0)
      _histJetAzimuth_pTmax180_->fill(dphi, event.weight());
  }


  // Finalize
  void D0_2004_S5992206::finalize() { 
    // Normalize histograms to unit area
    normalize(_histJetAzimuth_pTmax75_100);
    normalize(_histJetAzimuth_pTmax100_130);
    normalize(_histJetAzimuth_pTmax130_180);
    normalize(_histJetAzimuth_pTmax180_);
  }


}
