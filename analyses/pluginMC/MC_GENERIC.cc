// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Generic analysis looking at various distributions of final state particles
  class MC_GENERIC : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_GENERIC);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      const FinalState fs(Cuts::abseta < 5 && Cuts::pT > 500*MeV);
      declare(fs, "FS");
      declare(ChargedFinalState(fs), "CFS");

      // Histograms
      /// @todo Choose E/pT ranged based on input energies... can't do anything about kin. cuts, though
      book(_histMult   ,"Mult", 100, -0.5, 199.5);
      book(_histMultCh ,"MultCh", 100, -0.5, 199.5);

      book(_histPt   ,"Pt", 300, 0, 30);
      book(_histPtCh ,"PtCh", 300, 0, 30);

      book(_histE   ,"E", 100, 0, 200);
      book(_histECh ,"ECh", 100, 0, 200);

      book(_histEtaSumEt ,"EtaSumEt", 25, 0, 5);

      book(_histEta    ,"Eta", 50, -5, 5);
      book(_histEtaCh  ,"EtaCh", 50, -5, 5);
      _tmphistEtaPlus = Histo1D(25, 0, 5);
      _tmphistEtaMinus = Histo1D(25, 0, 5);
      _tmphistEtaChPlus = Histo1D(25, 0, 5);
      _tmphistEtaChMinus = Histo1D(25, 0, 5);

      book(_histRapidity    ,"Rapidity", 50, -5, 5);
      book(_histRapidityCh  ,"RapidityCh", 50, -5, 5);
      _tmphistRapPlus = Histo1D(25, 0, 5);
      _tmphistRapMinus = Histo1D(25, 0, 5);
      _tmphistRapChPlus = Histo1D(25, 0, 5);
      _tmphistRapChMinus = Histo1D(25, 0, 5);

      book(_histPhi    ,"Phi", 50, 0, TWOPI);
      book(_histPhiCh  ,"PhiCh", 50, 0, TWOPI);

      _histEtaPMRatio = bookScatter2D("EtaPMRatio");
      _histEtaChPMRatio = bookScatter2D("EtaChPMRatio");
      _histRapidityPMRatio = bookScatter2D("RapidityPMRatio");
      _histRapidityChPMRatio = bookScatter2D("RapidityChPMRatio");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      // Charged + neutral final state
      const FinalState& fs = apply<FinalState>(event, "FS");
      MSG_DEBUG("Total multiplicity = " << fs.size());
      _histMult->fill(fs.size(), weight);
      for (const Particle& p : fs.particles()) {
        _histEta->fill(p.eta(), weight);
        _histEtaSumEt->fill(p.abseta(), p.Et(), weight);
        (p.eta() > 0 ? _tmphistEtaPlus : _tmphistEtaMinus).fill(p.abseta(), weight);
        //
        _histRapidity->fill(p.rap(), weight);
        (p.rap() > 0 ? _tmphistRapPlus : _tmphistRapMinus).fill(p.absrap(), weight);
        //
        _histPt->fill(p.pT()/GeV, weight);
        _histE->fill(p.E()/GeV, weight);
        _histPhi->fill(p.phi(), weight);
      }

      // Same for the charged FS particles only
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histMultCh->fill(cfs.size(), weight);
      for (const Particle& p : cfs.particles()) {
        _histEtaCh->fill(p.eta(), weight);
        (p.eta() > 0 ? _tmphistEtaChPlus : _tmphistEtaChMinus).fill(p.abseta(), weight);
        //
        _histRapidityCh->fill(p.rap(), weight);
        (p.rap() > 0 ? _tmphistRapChPlus : _tmphistRapChMinus).fill(p.absrap(), weight);
        //
        _histPtCh->fill(p.pT()/GeV, weight);
        _histECh->fill(p.E()/GeV, weight);
        _histPhiCh->fill(p.phi(), weight);
      }

    }


    /// Finalize
    void finalize() {
      normalize({_histMult, _histEta, _histRapidity, _histPt, _histE, _histPhi});
      normalize({_histMultCh, _histEtaCh, _histRapidityCh, _histPtCh, _histECh, _histPhiCh});
      divide(_tmphistEtaPlus, _tmphistEtaMinus, _histEtaPMRatio);
      divide(_tmphistEtaChPlus, _tmphistEtaChMinus, _histEtaChPMRatio);
      divide(_tmphistRapPlus, _tmphistRapMinus, _histRapidityPMRatio);
      divide(_tmphistRapChPlus, _tmphistRapChMinus, _histRapidityChPMRatio);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histMult, _histEta, _histRapidity, _histPt, _histE, _histPhi;
    Histo1DPtr _histMultCh,  _histEtaCh, _histRapidityCh, _histPtCh, _histECh, _histPhiCh;
    Profile1DPtr _histEtaSumEt;
    Scatter2DPtr _histEtaPMRatio, _histEtaChPMRatio, _histRapidityPMRatio, _histRapidityChPMRatio;
    //@}

    /// @name Temporary histos used to calculate +/- rapidity ratio plots
    //@{
    Histo1D _tmphistEtaPlus, _tmphistEtaMinus, _tmphistEtaChPlus, _tmphistEtaChMinus;
    Histo1D _tmphistRapPlus, _tmphistRapMinus, _tmphistRapChPlus, _tmphistRapChMinus;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_GENERIC);

}
