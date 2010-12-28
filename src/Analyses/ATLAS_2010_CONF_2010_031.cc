// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {


  /// @brief ATLAS minimum bias analysis at 900/7000 GeV with nCh>=6
  class ATLAS_2010_CONF_2010_031 : public Analysis {
  public:

    ATLAS_2010_CONF_2010_031()
      : Analysis("ATLAS_2010_CONF_2010_031"),
        _Nevt_after_cuts(0.0)
    {
      setNeedsCrossSection(false);
    }


    void init() {
      ChargedFinalState cfs(-2.5, 2.5, 0.5*GeV);
      addProjection(cfs, "CFS");

//      if (sqrtS()/GeV >= 890 && sqrtS()/GeV <= 910) {
//        _h_dNch_deta = bookHistogram1D(2, 1, 1);
//        _h_dNevt_dNch = bookHistogram1D(2, 2, 1);
//        _h_dNch_dpT = bookHistogram1D(2, 3, 1);
//        _p_meanpT_Nch = bookProfile1D(2, 4, 1);
//      }
      if (sqrtS()/GeV >= 6990 && sqrtS()/GeV <= 7010) {
//        _h_dNch_deta = bookHistogram1D(3, 1, 1);
        _h_dNevt_dNch = bookHistogram1D(3, 2, 1);
//        _h_dNch_dpT = bookHistogram1D(3, 3, 1);
        _p_meanpT_Nch = bookProfile1D(3, 4, 1);
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      const ChargedFinalState& charged = applyProjection<ChargedFinalState>(event, "CFS");
      if (charged.size() < 1) {
        vetoEvent;
      }
      _Nevt_after_cuts += weight;

      _h_dNevt_dNch->fill(charged.size(), weight);
      foreach (const Particle& p, charged.particles()) {
        double pT = p.momentum().pT()/GeV;
//        _h_dNch_deta->fill(p.momentum().eta(), weight);
//        _h_dNch_dpT->fill(pT, weight/pT);
        _p_meanpT_Nch->fill(charged.size(), pT, weight);
      }
    }


    void finalize() {
//      double deta = 5.0;
//      scale(_h_dNch_deta, 1.0/_Nevt_after_cuts);
//      scale(_h_dNch_dpT, 1.0/_Nevt_after_cuts/TWOPI/deta);
      scale(_h_dNevt_dNch, 1.0/_Nevt_after_cuts);
    }


  private:

    AIDA::IHistogram1D* _h_dNch_deta;
    AIDA::IHistogram1D* _h_dNch_dpT;
    AIDA::IHistogram1D* _h_dNevt_dNch;
    AIDA::IProfile1D*   _p_meanpT_Nch;

    double _Nevt_after_cuts;

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<ATLAS_2010_CONF_2010_031> plugin_ATLAS_2010_CONF_2010_031;

}

