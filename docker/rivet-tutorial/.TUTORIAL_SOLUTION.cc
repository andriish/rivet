// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"

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

      // Projections
      FinalState leps(Cuts::pT > 20*GeV && Cuts::abspid == (_eemode ? PID::ELECTRON : PID::MUON));
      declare("Leps", leps);
      DressedLeptons dleps(0.1, Cuts::pT > 20*GeV && Cuts::abspid == (_eemode ? PID::ELECTRON : PID::MUON));
      declare("DLeps", dleps);
      FastJets jets(FinalState(), FastJets::ANTIKT, 0.4);
      declare("Jets", jets);

      // Book histograms
      vector<double> mll_bins = { 66., 74., 78., 82., 84., 86., 88., 89., 90., 91.,
                                  92., 93., 94., 96., 98., 100., 104., 108., 116. };
      book(_h["mll_bare"],    "mll_bare", mll_bins);
      book(_h["mll_dressed"], "mll_dressed", mll_bins);
      book(_h["njet_all"],    "njet_all",   6, -0.5,  5.5);
      book(_h["njet_iso"],    "njet_iso",   6, -0.5,  5.5);
      book(_h["nbjet_all"],   "nbjet_all",  3, -0.5,  2.5);
      book(_h["nbjet_iso"],   "nbjet_all",  3, -0.5,  2.5);
      book(_h["jet1pt"],      "jet1pt"      6,  20., 110.);
      book(_h["HT"],          "HT",         6,  20., 150.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // TASKS:
      // 1. Reconstruct the dilepton invariant mass to fill the histogram: difference between e/mu?
      // 2. Construct an alternative dilepton invariant mass histogram using dressed leptons
      // 3. Find jets -- what muon treatment? -- and histogram multiplicity and HT: e/mu difference?
      // 4. Isolate the jets from leptons and look again
      // 5. Count the number of b-jets, using the ATLAS pT_B > 5 GeV tag definition

      // Get leptons (bare & dressed)
      const Particles leps = apply<FinalState>(event, "Leps").particlesByPt();
      const Particles dleps = apply<FinalState>(event, "DLeps").particlesByPt();

      // Get jets and remove overlaps
      const Jets alljets = apply<JetFinder>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::absrap < 4.5);
      const Jets isojets = discardIfAnyDeltaRLess(alljets, dleps, 0.4);

      // m_ll
      if (leps.size() == 2) _h["mll_bare"]->fill( (leps[0].mom() + leps[1].mom()).mass()/GeV );
      if (dleps.size() == 2) _h["mll_dressed"]->fill( (dleps[0].mom() + dleps[1].mom()).mass()/GeV );

      // Jets
      _h["njet_all"]->fill(alljets.size());
      _h["njet_iso"]->fill(isojets.size());
      _h["jet1pt"]->fill(isojets[0].pT()/GeV);
      _h["HT"]->fill(sum(isojets, pT, 0.0)/GeV + sum(head(dleps, 2), pT, 0.0)/GeV);

      // b-jets
      _h["nbjet_all"]->fill(select(alljets, hasBTag(Cuts::pT > 5*GeV)).size());
      _h["nbjet_iso"]->fill(select(isojets, hasBTag(Cuts::pT > 5*GeV)).size());

    }


    /// Post-run histogram manipulations
    void finalize() {

      // Normalise histograms
      const double sf = crossSection()/picobarn / sumW();
      for (auto hist : _h) { scale(hist.second, sf); }
    }


    /// Histograms
    map<string, Histo1DPtr> _h;

    /// Mode flag
    bool _eemode;

  };


  DECLARE_RIVET_PLUGIN(TUTORIAL);

}
