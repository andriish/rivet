// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetYODA.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z + jets events
  class MC_ZJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_ZJETS()
      : MC_JetAnalysis("MC_ZJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      ZFinder zfinder(fs, -3.5, 3.5, 25.0*GeV, ELECTRON, 65.0*GeV, 115.0*GeV, 0.2, true, true);
      addProjection(zfinder, "ZFinder");
      FastJets jetpro(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.4);
      addProjection(jetpro, "Jets");

      _h_Z_jet1_deta = bookHisto1D("Z_jet1_deta", 50, -5.0, 5.0);
      _h_Z_jet1_dR = bookHisto1D("Z_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size()!=1) {
        vetoEvent;
      }
      const double weight = e.weight();

      FourMomentum zmom(zfinder.bosons()[0].momentum());

      const Jets& jets = applyProjection<FastJets>(e, "Jets").jetsByPt(m_jetptcut);
      if (jets.size() > 0) {
        _h_Z_jet1_deta->fill(zmom.eta()-jets[0].momentum().eta(), weight);
        _h_Z_jet1_dR->fill(deltaR(zmom, jets[0].momentum()), weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      scale(_h_Z_jet1_deta, crossSection()/sumOfWeights());
      scale(_h_Z_jet1_dR, crossSection()/sumOfWeights());

      MC_JetAnalysis::finalize();
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Z_jet1_deta;
    Histo1DPtr _h_Z_jet1_dR;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ZJETS);

}
