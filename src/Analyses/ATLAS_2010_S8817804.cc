// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief ATLAS inclusive jet pT spectrum, di-jet Mass and di-jet chi
  class ATLAS_2010_S8817804: public Analysis {
  public:

    ATLAS_2010_S8817804() : Analysis("ATLAS_2010_S8817804")
    {
      setBeams(PROTON, PROTON);
      setNeedsCrossSection(true);
      _sumW = 0;
    }


  private:

    enum Alg { AKT4=0, AKT6=1 };


  public:

    void init() {
      FinalState fs;
      addProjection(fs, "FinalState");
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.6), "AntiKT06");
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.4), "AntiKT04");

      double ybins[] = { 0.0, 0.3, 0.8, 1.2, 2.1, 2.8 };
      double massBinsForChi[] = { 340, 520, 800, 1200 };

      for (size_t alg = 0; alg < 2; ++alg) {
        int ptDsOffset = 1;
        for (size_t i = 0; i < 5; ++i) {
          _pTHistos[alg].addHistogram(ybins[i], ybins[i+1], bookHistogram1D(i + ptDsOffset, 1, 1));
        }
        ptDsOffset += 5;

        int massDsOffset = 11;
        for (size_t i = 0; i < 5; ++i) {
          _massVsY[alg].addHistogram(ybins[i], ybins[i+1], bookHistogram1D(i + massDsOffset, 1, 1));
        }
        massDsOffset += 5;

        int chiDsOffset = 21;
        for (size_t i = 0; i < 3; ++i) {
          _chiVsMass[alg].addHistogram(massBinsForChi[i], massBinsForChi[i+1], bookHistogram1D(i + chiDsOffset, 1, 1));
        }
        chiDsOffset += 3;
      }
    }


    void analyze(const Event& evt) {
      Jets jetAr[2];
      jetAr[AKT4] = applyProjection<FastJets>(evt, "AntiKT06").jetsByPt(30*GeV);
      jetAr[AKT6] = applyProjection<FastJets>(evt, "AntiKT04").jetsByPt(30*GeV);

      cout << -1 << endl;

      // Identify the dijets
      for (size_t alg = 0; alg < 2; ++alg) {
        cout << 0 << endl;

        vector<FourMomentum> leadjets;
        foreach (const Jet& jet, jetAr[alg]) {
          cout << 0.1 << endl;
          double pT = jet.momentum().pT();
          double absy = fabs(jet.momentum().rapidity());
          cout << 0.2 << endl;
          _pTHistos[alg].fill(absy, pT/GeV, evt.weight());

          cout << 0.3 << endl;

          if (absy < 2.8 && leadjets.size() < 2) {
            if (leadjets.empty() && pT < 60*GeV) continue;
            leadjets.push_back(jet.momentum());
          }
          cout << 0.4 << endl;
        }

        cout << 1 << endl;

        // Veto if no acceptable dijet found
        if (leadjets.size() < 2) {
          MSG_DEBUG("Could not find two suitable leading jets");
          continue;
        }
        _sumW += evt.weight();

        cout << 2 << endl;

        const double rap1 = leadjets[0].rapidity();
        const double rap2 = leadjets[1].rapidity();
        const double mass = (leadjets[0] + leadjets[1]).mass();
        //const double ymax = (fabs(rap1) > fabs(rap2)) ? fabs(rap1) : fabs(rap2);
        const double ymax = max(fabs(rap1), fabs(rap2));
        const double chi = exp(fabs(rap1 - rap2));
        if (fabs(rap1 + rap2) < 2.2) {
          _chiVsMass[alg].fill(mass/GeV, chi, evt.weight());
        }
        _massVsY[alg].fill(ymax, mass/GeV, evt.weight());

      }
    }


    void finalize() {
      for (size_t alg = 0; alg < 2; ++alg) {
        _pTHistos[alg].scale(crossSectionPerEvent()/picobarn, this);
        _massVsY[alg].scale(crossSection()/_sumW/picobarn, this);
        _chiVsMass[alg].scale(crossSection()/_sumW/picobarn, this);
      }
    }


  private:

    /// Counter for weights passing the dijet requirement cut
    double _sumW;

    /// The inclusive pT spectrum for akt6 and akt4 jets (array index is jet type from enum above)
    BinnedHistogram<double> _pTHistos[2];

    /// The di-jet mass spectrum binned in rapidity for akt6 and akt4 jets (array index is jet type from enum above)
    BinnedHistogram<double> _massVsY[2];

    /// The di-jet chi distribution binned in mass for akt6 and akt4 jets (array index is jet type from enum above)
    BinnedHistogram<double> _chiVsMass[2];

  };


  AnalysisBuilder<ATLAS_2010_S8817804> plugin_ATLAS_2010_S8817804;

}
