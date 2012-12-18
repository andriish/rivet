// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// TOTEM elastic and total cross-section measurement
  class TOTEM_2012_002 : public Analysis {
  public:

    TOTEM_2012_002()
      : Analysis("TOTEM_2012_002")
    {    }


  public:

    void init() {
      addProjection(ChargedFinalState(), "CFS");
      _hist_tlow  = bookHistogram1D(1, 1, 1);
      _hist_thigh = bookHistogram1D(2, 1, 1);
      _hist_sigma = bookHistogram1D(3, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
      if (cfs.size() > 2) MSG_WARNING("Final state includes more than two charged particles !");
      _hist_sigma->fill(sqrtS()/GeV, weight);

      foreach (const Particle& p, cfs.particles()) {
        if (p.momentum().eta() > 0. && p.pdgId() == PROTON) {
          double t = sqr(p.momentum().pT());
          _hist_tlow->fill(t, weight);
          _hist_thigh->fill(t, weight);
        }
      }
    }


    void finalize() {
      scale(_hist_tlow, crossSection()/millibarn/sumOfWeights());
      scale(_hist_thigh, crossSection()/millibarn/sumOfWeights());
      scale(_hist_sigma, crossSection()/millibarn/sumOfWeights());
    }


  private:

    AIDA::IHistogram1D *_hist_tlow;
    AIDA::IHistogram1D *_hist_thigh;
    AIDA::IHistogram1D *_hist_sigma;

  };


  DECLARE_RIVET_PLUGIN(TOTEM_2012_002);

}
