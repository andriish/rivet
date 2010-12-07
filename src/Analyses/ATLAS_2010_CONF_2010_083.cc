// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "LWH/Histogram1D.h"

namespace Rivet {

  class ATLAS_2010_CONF_2010_083 : public Analysis {
  public:

    ATLAS_2010_CONF_2010_083() : Analysis("ATLAS_2010_CONF_2010_083") {
      setNeedsCrossSection(true);
    }


    void init() {
      FinalState fs(-4.5, 4.5, 0.0*GeV);
      FastJets jetsproj(fs, FastJets::ANTIKT, 0.6);
      addProjection(jetsproj, "Jets");

      for (size_t i=0 ; i<4 ; i++) {
        _h_dphi[i] = bookHistogram1D(i+1, 1, 1);
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      const FastJets & jetsproj = applyProjection<FastJets>(event, "Jets");
      Jets alljets;
      foreach (const Jet jet, jetsproj.jetsByPt(100.0*GeV)) {
        if (fabs(jet.momentum().rapidity())<2.8) {
          alljets.push_back(jet);
        }
      }
      Jets jets;
      foreach (const Jet jet, alljets) {
        if (fabs(jet.momentum().rapidity())<0.8) {
          jets.push_back(jet);
        }
        if (jets.size()==2) break;
      }

      if (jets.size() < 2 || jets[0].momentum().pT() < 110.0*GeV) {
        vetoEvent;
      }

      const double dphi = deltaPhi(jets[0].momentum().phi(), jets[1].momentum().phi());
      const double ptlead = jets[0].momentum().pT();

      if (ptlead > 310.0*GeV) {
        _h_dphi[3]->fill(dphi, weight);
      }
      else if (ptlead > 210.0*GeV) {
        _h_dphi[2]->fill(dphi, weight);
      }
      else if (ptlead > 160.0*GeV) {
        _h_dphi[1]->fill(dphi, weight);
      }
      else {
        _h_dphi[0]->fill(dphi, weight);
      }
    }

    void finalize() {
      for (size_t i=0 ; i<4 ; i++) {
        normalize(_h_dphi[i], 1.0);
      }
    }


  private:

    AIDA::IHistogram1D* _h_dphi[4];
  };

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<ATLAS_2010_CONF_2010_083> plugin_ATLAS_2010_CONF_2010_083;

}

