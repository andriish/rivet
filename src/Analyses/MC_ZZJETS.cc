// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/RivetYODA.hh"

namespace Rivet {

  /// @brief MC validation analysis for Z[ee]Z[mumu] + jets events
  class MC_ZZJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_ZZJETS()
      : MC_JetAnalysis("MC_ZZJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      ZFinder zeefinder(fs, -3.5, 3.5, 25.0*GeV, PID::ELECTRON, 65.0*GeV, 115.0*GeV, 0.2, true, true);
      addProjection(zeefinder, "ZeeFinder");

      VetoedFinalState zmminput;
      zmminput.addVetoOnThisFinalState(zeefinder);
      ZFinder zmmfinder(zmminput, -3.5, 3.5, 25.0*GeV, PID::MUON, 65.0*GeV, 115.0*GeV, 0.2, true, true);
      addProjection(zmmfinder, "ZmmFinder");

      VetoedFinalState jetinput;
      jetinput
          .addVetoOnThisFinalState(zeefinder)
          .addVetoOnThisFinalState(zmmfinder);
      FastJets jetpro(jetinput, FastJets::ANTIKT, 0.4);
      addProjection(jetpro, "Jets");

      // correlations with jets
      _h_ZZ_jet1_deta = bookHisto1D("ZZ_jet1_deta", 70, -7.0, 7.0);
      _h_ZZ_jet1_dR = bookHisto1D("ZZ_jet1_dR", 25, 1.5, 7.0);
      _h_Ze_jet1_dR = bookHisto1D("Ze_jet1_dR", 25, 0.0, 7.0);

      // global stuff
      _h_HT = bookHisto1D("HT", logspace(100, 100.0, 0.5*sqrtS()));

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = e.weight();

      const ZFinder& zeefinder = applyProjection<ZFinder>(e, "ZeeFinder");
      if (zeefinder.bosons().size()!=1) {
        vetoEvent;
      }

      const ZFinder& zmmfinder = applyProjection<ZFinder>(e, "ZmmFinder");
      if (zmmfinder.bosons().size()!=1) {
        vetoEvent;
      }

      FourMomentum zee(zeefinder.bosons()[0].momentum());
      FourMomentum zmm(zmmfinder.bosons()[0].momentum());
      FourMomentum zz(zee+zmm);
      // find leptons
      FourMomentum ep(zeefinder.constituents()[0].momentum()),
        em(zeefinder.constituents()[1].momentum()),
        mp(zmmfinder.constituents()[0].momentum()),
        mm(zmmfinder.constituents()[1].momentum());

      const Jets& jets = applyProjection<FastJets>(e, "Jets").jetsByPt(m_jetptcut);
      if (jets.size() > 0) {
        _h_ZZ_jet1_deta->fill(zz.eta()-jets[0].momentum().eta(), weight);
        _h_ZZ_jet1_dR->fill(deltaR(zz, jets[0].momentum()), weight);
        _h_Ze_jet1_dR->fill(deltaR(ep, jets[0].momentum()), weight);
      }

      double HT=ep.pT()+em.pT()+mp.pT()+mm.pT();
      foreach (const Jet& jet, jets) {
        HT+=jet.momentum().pT();
      }
      if (HT>0.0) _h_HT->fill(HT, weight);

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      double norm=crossSection()/sumOfWeights();
      scale(_h_ZZ_jet1_deta, norm);
      scale(_h_ZZ_jet1_dR, norm);
      scale(_h_Ze_jet1_dR, norm);
      scale(_h_HT, norm);

      MC_JetAnalysis::finalize();
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_ZZ_jet1_deta;
    Histo1DPtr _h_ZZ_jet1_dR;
    Histo1DPtr _h_Ze_jet1_dR;
    Histo1DPtr _h_HT;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ZZJETS);

}
