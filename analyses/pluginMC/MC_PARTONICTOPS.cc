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
      declare(PartonicTops(PartonicTops::DecayMode::ALL), "AllTops");
      declare(PartonicTops(PartonicTops::DecayMode::ALL, true, false, Cuts::OPEN, PartonicTops::WhichTop::FIRST), "AllTopsFirst"); ///< @todo API ick!
      declare(PartonicTops(PartonicTops::DecayMode::E_MU), "LeptonicTops");
      declare(PartonicTops(PartonicTops::DecayMode::HADRONIC), "HadronicTops");

      // Book histograms
      book(_h_tall_n,  "t_all_n", linspace(5, -0.5, 4.5));
      book(_h_tall_pt, "t_all_pT", logspace(50, 1, 500));
      book(_h_tall_y,  "t_all_y", linspace(50, -5, 5));

      book(_h_tall_n_first,  "t_all_n_firsttop", linspace(5, -0.5, 4.5));
      book(_h_tall_pt_first, "t_all_pT_firsttop", logspace(50, 1, 500));
      book(_h_tall_y_first,  "t_all_y_firsttop", linspace(50, -5, 5));

      book(_h_tlep_n,  "t_lep_n", linspace(5, -0.5, 4.5));
      book(_h_tlep_pt, "t_lep_pT", logspace(50, 1, 500));
      book(_h_tlep_y,  "t_lep_y", linspace(50, -5, 5));

      book(_h_thad_n,  "t_had_n", linspace(5, -0.5, 4.5));
      book(_h_thad_pt, "t_had_pT", logspace(50, 1, 500));
      book(_h_thad_y,  "t_had_y", linspace(50, -5, 5));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& alltops = apply<PartonicTops>(event, "AllTops").particlesByPt();
      _h_tall_n->fill(alltops.size());
      for (const Particle& t : alltops) {
        _h_tall_pt->fill(t.pT()/GeV);
        _h_tall_y->fill(t.rap());
      }

      const Particles& alltops_first = apply<PartonicTops>(event, "AllTopsFirst").particlesByPt();
      _h_tall_n_first->fill(alltops_first.size());
      for (const Particle& t : alltops_first) {
        _h_tall_pt_first->fill(t.pT()/GeV);
        _h_tall_y_first->fill(t.rap());
      }

      const Particles& leptops = apply<PartonicTops>(event, "LeptonicTops").particlesByPt();
      _h_tlep_n->fill(leptops.size());
      for (const Particle& t : leptops) {
        _h_tlep_pt->fill(t.pT()/GeV);
        _h_tlep_y->fill(t.rap());
      }

      const Particles& hadtops = apply<PartonicTops>(event, "HadronicTops").particlesByPt();
      _h_thad_n->fill(hadtops.size());
      for (const Particle& t : hadtops) {
        _h_thad_pt->fill(t.pT()/GeV);
        _h_thad_y->fill(t.rap());
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize({_h_tall_n,  _h_tall_n_first, _h_tlep_n, _h_thad_n});
      normalize({_h_tall_pt, _h_tall_pt_first, _h_tlep_pt, _h_thad_pt});
      normalize({_h_tall_y,  _h_tall_y_first, _h_tlep_y, _h_thad_y});
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_tall_n, _h_tall_n_first, _h_tlep_n, _h_thad_n;
    Histo1DPtr _h_tall_pt, _h_tall_pt_first, _h_tlep_pt, _h_thad_pt;
    Histo1DPtr _h_tall_y, _h_tall_y_first, _h_tlep_y, _h_thad_y;
    //@}


  };


  DECLARE_RIVET_PLUGIN(MC_PARTONICTOPS);

}
