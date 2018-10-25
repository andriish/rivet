// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  /// Find and plot partonic top properties (requires tops in event record)
  class MC_PARTONICTOPS : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_PARTONICTOPS);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(PartonicTops(PartonicTops::ALL), "AllTops");
      declare(PartonicTops(PartonicTops::E_MU), "LeptonicTops");
      declare(PartonicTops(PartonicTops::HADRONIC), "HadronicTops");

      // Book histograms
      _h_tall_n  = bookHisto1D("t_all_n", linspace(5, -0.5, 4.5));
      _h_tall_pt = bookHisto1D("t_all_pT", logspace(50, 1, 500));
      _h_tall_y  = bookHisto1D("t_all_y", linspace(50, -5, 5));

      _h_tlep_n  = bookHisto1D("t_lep_n", linspace(5, -0.5, 4.5));
      _h_tlep_pt = bookHisto1D("t_lep_pT", logspace(50, 1, 500));
      _h_tlep_y  = bookHisto1D("t_lep_y", linspace(50, -5, 5));

      _h_thad_n  = bookHisto1D("t_had_n", linspace(5, -0.5, 4.5));
      _h_thad_pt = bookHisto1D("t_had_pT", logspace(50, 1, 500));
      _h_thad_y  = bookHisto1D("t_had_y", linspace(50, -5, 5));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& alltops = apply<PartonicTops>(event, "AllTops").particlesByPt();
      _h_tall_n->fill(alltops.size(), event.weight());
      for (const Particle& t : alltops) {
        _h_tall_pt->fill(t.pT()/GeV, event.weight());
        _h_tall_y->fill(t.rap(), event.weight());
      }

      const Particles& leptops = apply<PartonicTops>(event, "LeptonicTops").particlesByPt();
      _h_tlep_n->fill(leptops.size(), event.weight());
      for (const Particle& t : leptops) {
        _h_tlep_pt->fill(t.pT()/GeV, event.weight());
        _h_tlep_y->fill(t.rap(), event.weight());
      }

      const Particles& hadtops = apply<PartonicTops>(event, "HadronicTops").particlesByPt();
      _h_thad_n->fill(hadtops.size(), event.weight());
      for (const Particle& t : hadtops) {
        _h_thad_pt->fill(t.pT()/GeV, event.weight());
        _h_thad_y->fill(t.rap(), event.weight());
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize({_h_tall_n, _h_tlep_n, _h_thad_n});
      normalize({_h_tall_pt, _h_tlep_pt, _h_thad_pt});
      normalize({_h_tall_y, _h_tlep_y, _h_thad_y});
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_tall_n, _h_tlep_n, _h_thad_n;
    Histo1DPtr _h_tall_pt, _h_tlep_pt, _h_thad_pt;
    Histo1DPtr _h_tall_y, _h_tlep_y, _h_thad_y;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_PARTONICTOPS);


}
