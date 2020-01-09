// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Tutorial analysis with dilepton and jet observables
  class TUTORIAL : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TUTORIAL);


    /// Book histograms and initialise projections before the run
    void init() {

      // Electron or muon mode? (default to electron)
      _eemode = (getOption("LMODE") != "MU");

      // Book histograms
      vector<double> mll_bins = { 66., 74., 78., 82., 84., 86., 88., 89., 90., 91.,
                                  92., 93., 94., 96., 98., 100., 104., 108., 116. };
      book(_h["mll"], "mass_ll", mll_bins);
      //book(_h["jets_excl"],  "jets_excl",   6, -0.5,  5.5);
      //book(_h["bjets_excl"], "bjets_excl",  3, -0.5,  2.5);
      //book(_h["HT"],         "HT",          6,  20., 110.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // TASKS:
      // 1. Reconstruct the dilepton invariant mass to fill the histogram: difference between e/mu?
      // 2. Construct an alternative dilepton invariant mass histogram using dressed leptons
      // 3. Find jets -- what muon treatment? -- and histogram multiplicity and HT: e/mu difference?
      // 4. Isolate the jets from leptons and look again
      // 5. Count the number of b-jets, using the ATLAS pT_B > 5 GeV tag definition

      // ...
      _h["mll"]->fill(mll);
    }


    /// Post-run histogram manipulations
    void finalize() {

      // Normalise histograms
      const double sf = crossSection() / sumW();
      for (auto hist : _h) { scale(hist.second, sf); }
    }


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    bool _eemode;
    //@}

  };


  DECLARE_RIVET_PLUGIN(TUTORIAL);

}
