// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  /// @brief MC validation analysis for W + jets events
  class MC_WJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_WJETS()
      : MC_JetAnalysis("MC_WJETS", 4, "Jets")
    {
      setNeedsCrossSection(true);
    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      WFinder wfinder(-3.5, 3.5, 25.0*GeV, ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      addProjection(wfinder, "WFinder");
      FastJets jetpro(wfinder.remainingFinalState(), FastJets::KT, 0.7);
      addProjection(jetpro, "Jets");

      _h_W_mass = bookHistogram1D("W_mass", 50, 55.0, 105.0);
      _h_W_pT = bookHistogram1D("W_pT", logBinEdges(100, 1.0, 0.5*sqrtS()));
      _h_W_pT_peak = bookHistogram1D("W_pT_peak", 25, 0.0, 25.0);
      _h_W_y = bookHistogram1D("W_y", 40, -4.0, 4.0);
      _h_W_phi = bookHistogram1D("W_phi", 25, 0.0, TWOPI);
      _h_W_jet1_deta = bookHistogram1D("W_jet1_deta", 50, -5.0, 5.0);
      _h_W_jet1_dR = bookHistogram1D("W_jet1_dR", 25, 0.5, 7.0);
      _h_lepton_pT = bookHistogram1D("lepton_pT", logBinEdges(100, 10.0, 0.25*sqrtS()));
      _h_lepton_eta = bookHistogram1D("lepton_eta", 40, -4.0, 4.0);

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const WFinder& wfinder = applyProjection<WFinder>(e, "WFinder");
      if (wfinder.particles().size()!=1) {
        vetoEvent;
      }
      const double weight = e.weight();

      FourMomentum wmom(wfinder.particles()[0].momentum());
      _h_W_mass->fill(wmom.mass(),weight);
      _h_W_pT->fill(wmom.pT(),weight);
      _h_W_pT_peak->fill(wmom.pT(),weight);
      _h_W_y->fill(wmom.rapidity(),weight);
      _h_W_phi->fill(wmom.azimuthalAngle(),weight);
      foreach (const Particle& l, wfinder.constituentLeptonsFinalState().particles()) {
        _h_lepton_pT->fill(l.momentum().pT(), weight);
        _h_lepton_eta->fill(l.momentum().eta(), weight);
      }

      const FastJets& jetpro = applyProjection<FastJets>(e, "Jets");
      const Jets& jets = jetpro.jetsByPt(20.0*GeV);
      if (jets.size() > 0) {
        _h_W_jet1_deta->fill(wmom.eta()-jets[0].momentum().eta(), weight);
        _h_W_jet1_dR->fill(deltaR(wmom, jets[0].momentum()), weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      scale(_h_W_mass, crossSection()/sumOfWeights());
      scale(_h_W_pT, crossSection()/sumOfWeights());
      scale(_h_W_pT_peak, crossSection()/sumOfWeights());
      scale(_h_W_y, crossSection()/sumOfWeights());
      scale(_h_W_phi, crossSection()/sumOfWeights());
      scale(_h_W_jet1_deta, crossSection()/sumOfWeights());
      scale(_h_W_jet1_dR, crossSection()/sumOfWeights());
      scale(_h_lepton_pT, crossSection()/sumOfWeights());
      scale(_h_lepton_eta, crossSection()/sumOfWeights());

      MC_JetAnalysis::finalize();
    }

    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_W_mass;
    AIDA::IHistogram1D * _h_W_pT;
    AIDA::IHistogram1D * _h_W_pT_peak;
    AIDA::IHistogram1D * _h_W_y;
    AIDA::IHistogram1D * _h_W_phi;
    AIDA::IHistogram1D * _h_W_jet1_deta;
    AIDA::IHistogram1D * _h_W_jet1_dR;
    AIDA::IHistogram1D * _h_lepton_pT;
    AIDA::IHistogram1D * _h_lepton_eta;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_WJETS> plugin_MC_WJETS;

}
