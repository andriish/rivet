// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/D0_2009_S8349509.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  D0_2009_S8349509::D0_2009_S8349509() : Analysis("D0_2009_S8349509") {
    setBeams(PROTON, ANTIPROTON);

    ZFinder zfinder(-1.7, 1.7, 15.0*GeV, MUON, 65.0*GeV, 115.0*GeV, 0.2);
    addProjection(zfinder, "ZFinder");

    FastJets conefinder(zfinder.remainingFinalState(), FastJets::D0ILCONE, 0.5, 20.0*GeV);
    addProjection(conefinder, "ConeFinder");
  }


  void D0_2009_S8349509::init() {

    _h_dphi_jet_Z25 = bookHistogram1D(1, 1, 1);
    _h_dphi_jet_Z45 = bookHistogram1D(2, 1, 1);

    _h_dy_jet_Z25 = bookHistogram1D(3, 1, 1);
    _h_dy_jet_Z45 = bookHistogram1D(4, 1, 1);

    _h_yboost_jet_Z25 = bookHistogram1D(5, 1, 1);
    _h_yboost_jet_Z45 = bookHistogram1D(6, 1, 1);
    
    _inclusive_Z_sumofweights = 0.0;
  }


  void D0_2009_S8349509::analyze(const Event& event) {

    double weight = event.weight();

    const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");
    if (zfinder.particles().size()==1) {
      // count inclusive sum of weights for histogram normalisation
      _inclusive_Z_sumofweights += weight;
      
      Jets jets;
      foreach (const Jet& j, applyProjection<JetAlg>(event, "ConeFinder").jetsByPt()) {
        if (fabs(j.momentum().pseudorapidity()) < 2.8) {
          jets.push_back(j);
          break;
        }
      }
      
      // Return if there are no jets:
      if(jets.size()<1) {
        getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
            << " because no jets pass cuts " << endl;
        vetoEvent;
      }
      
      // cut on Delta R between jet and muons
      foreach (const Jet& j, jets) {
        foreach (const Particle& mu, zfinder.constituentsFinalState().particles()) {
          if (deltaR(mu.momentum(), j.momentum()) < 0.5) {
            vetoEvent;
          }
        }
      }
      
      const FourMomentum Zmom = zfinder.particles()[0].momentum();
      const FourMomentum jetmom = jets[0].momentum();
      double yZ = Zmom.rapidity();
      double yjet = jetmom.rapidity();
      double dphi = deltaPhi(Zmom.phi(), jetmom.phi());
      double dy = fabs(yZ-yjet);
      double yboost = fabs(yZ+yjet)/2.0;

      if (Zmom.pT()>25.0*GeV) {
        _h_dphi_jet_Z25->fill(dphi,weight);
        _h_dy_jet_Z25->fill(dy, weight);
        _h_yboost_jet_Z25->fill(yboost, weight);
      }
      if (Zmom.pT()>45.0*GeV) {
        _h_dphi_jet_Z45->fill(dphi,weight);
        _h_dy_jet_Z45->fill(dy, weight);
        _h_yboost_jet_Z45->fill(yboost, weight);
      }
    }

  }


  void D0_2009_S8349509::finalize() {
    if (_inclusive_Z_sumofweights==0.0) return;
    scale(_h_dphi_jet_Z25, 1.0/_inclusive_Z_sumofweights);
    scale(_h_dphi_jet_Z45, 1.0/_inclusive_Z_sumofweights);
    scale(_h_dy_jet_Z25, 1.0/_inclusive_Z_sumofweights);
    scale(_h_dy_jet_Z45, 1.0/_inclusive_Z_sumofweights);
    scale(_h_yboost_jet_Z25, 1.0/_inclusive_Z_sumofweights);
    scale(_h_yboost_jet_Z45, 1.0/_inclusive_Z_sumofweights);
  }



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<D0_2009_S8349509> plugin_D0_2009_S8349509;

}
