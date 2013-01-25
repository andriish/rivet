// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  /// @brief MC validation analysis for W^+[enu]W^-[munu] + jets events
  class MC_WWJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_WWJETS()
      : MC_JetAnalysis("MC_WWJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      WFinder wenufinder(fs, -3.5, 3.5, 25.0*GeV, ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      addProjection(wenufinder, "WenuFinder");

      VetoedFinalState wmnuinput;
      wmnuinput.addVetoOnThisFinalState(wenufinder);
      WFinder wmnufinder(wmnuinput, -3.5, 3.5, 25.0*GeV, MUON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      addProjection(wmnufinder, "WmnuFinder");

      VetoedFinalState jetinput;
      jetinput
          .addVetoOnThisFinalState(wenufinder)
          .addVetoOnThisFinalState(wmnufinder);
      FastJets jetpro(jetinput, FastJets::ANTIKT, 0.4);
      addProjection(jetpro, "Jets");

      // correlations with jets
      _h_WW_jet1_deta = bookHistogram1D("WW_jet1_deta", 70, -7.0, 7.0);
      _h_WW_jet1_dR = bookHistogram1D("WW_jet1_dR", 25, 1.5, 7.0);
      _h_We_jet1_dR = bookHistogram1D("We_jet1_dR", 25, 0.0, 7.0);

      // global stuff
      _h_HT = bookHistogram1D("HT", logBinEdges(100, 100.0, 0.5*sqrtS()));
      _h_jets_m_12 = bookHistogram1D("jets_m_12", logBinEdges(100, 1.0, 0.25*sqrtS()));

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = e.weight();

      const WFinder& wenufinder = applyProjection<WFinder>(e, "WenuFinder");
      if (wenufinder.bosons().size()!=1) {
        vetoEvent;
      }

      const WFinder& wmnufinder = applyProjection<WFinder>(e, "WmnuFinder");
      if (wmnufinder.bosons().size()!=1) {
        vetoEvent;
      }

      FourMomentum wenu(wenufinder.bosons()[0].momentum());
      FourMomentum wmnu(wmnufinder.bosons()[0].momentum());
      FourMomentum ww(wenu+wmnu);
      // find leptons
      FourMomentum ep=wenufinder.constituentLeptons()[0].momentum();
      FourMomentum enu=wenufinder.constituentNeutrinos()[0].momentum();
      FourMomentum mm=wmnufinder.constituentLeptons()[0].momentum();
      FourMomentum mnu=wmnufinder.constituentNeutrinos()[0].momentum();

      const Jets& jets = applyProjection<FastJets>(e, "Jets").jetsByPt(m_jetptcut);
      if (jets.size() > 0) {
        _h_WW_jet1_deta->fill(ww.eta()-jets[0].momentum().eta(), weight);
        _h_WW_jet1_dR->fill(deltaR(ww, jets[0].momentum()), weight);
        _h_We_jet1_dR->fill(deltaR(ep, jets[0].momentum()), weight);
      }

      double HT=ep.pT()+mm.pT()+FourMomentum(enu+mnu).pT();
      foreach (const Jet& jet, jets) {
        HT+=jet.momentum().pT();
      }
      if (HT>0.0) _h_HT->fill(HT, weight);

      if (jets.size()>1) {
        FourMomentum jet1(jets[0].momentum());
        FourMomentum jet2(jets[1].momentum());
        _h_jets_m_12->fill(FourMomentum(jet1+jet2).mass(), weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      double norm=crossSection()/sumOfWeights();
      scale(_h_WW_jet1_deta, norm);
      scale(_h_WW_jet1_dR, norm);
      scale(_h_We_jet1_dR, norm);
      scale(_h_jets_m_12, norm);
      scale(_h_HT, norm);

      MC_JetAnalysis::finalize();
    }

    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_WW_jet1_deta;
    AIDA::IHistogram1D * _h_WW_jet1_dR;
    AIDA::IHistogram1D * _h_We_jet1_dR;
    AIDA::IHistogram1D * _h_jets_m_12;
    AIDA::IHistogram1D * _h_HT;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_WWJETS);

}
