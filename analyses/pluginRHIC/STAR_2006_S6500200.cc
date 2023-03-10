// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// STAR identified hadron spectra in pp at 200 GeV
  class STAR_2006_S6500200 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2006_S6500200);


    /// Book projections and histograms
    void init() {
      ChargedFinalState bbc1(Cuts::etaIn(-5.0,-3.3)); // beam-beam-counter trigger
      ChargedFinalState bbc2(Cuts::etaIn( 3.3, 5.0)); // beam-beam-counter trigger
      declare(bbc1, "BBC1");
      declare(bbc2, "BBC2");

      IdentifiedFinalState pionfs(Cuts::abseta < 2.5 && Cuts::pT > 0.3*GeV);
      IdentifiedFinalState protonfs(Cuts::abseta < 2.5 && Cuts::pT > 0.4*GeV);
      pionfs.acceptIdPair(PID::PIPLUS);
      protonfs.acceptIdPair(PID::PROTON);
      declare(pionfs, "PionFS");
      declare(protonfs, "ProtonFS");

      book(_h_pT_piplus     ,2, 1, 1); // full range pion binning
      book(_h_pT_piminus    ,7, 1, 1); // full range pion binning
      book(_tmp_pT_piplus   ,"TMP/pT_piplus", refData(25, 1, 2)); // pi histo compatible with more restricted proton binning
      book(_tmp_pT_piminus  ,"TMP/pT_piminus", refData(26, 1, 2)); // pi histo compatible with more restricted proton binning
      book(_h_pT_proton     ,12, 1, 1);
      book(_h_pT_antiproton ,17, 1, 1);

      book(_s_piminus_piplus, 23, 1, 2);
      book(_s_antipr_pr     , 24, 1, 2);
      book(_s_pr_piplus     , 25, 1, 2);
      book(_s_antipr_piminus, 26, 1, 2);

      book(_sumWeightSelected, "_sumWeightSelected");
    }


    /// Do the analysis
    void analyze(const Event& event) {
      const ChargedFinalState& bbc1 = apply<ChargedFinalState>(event, "BBC1");
      const ChargedFinalState& bbc2 = apply<ChargedFinalState>(event, "BBC2");
      if (bbc1.size() < 1 || bbc2.size() < 1) {
        MSG_DEBUG("Failed beam-beam-counter trigger");
        vetoEvent;
      }

      const IdentifiedFinalState& pionfs = apply<IdentifiedFinalState>(event, "PionFS");
      for (const Particle& p : pionfs.particles()) {
        if (p.absrap() < 0.5) {
          /// @todo Use a binned counter to avoid this bin width cancellation hack
          const double pT = p.pT() / GeV;
          ((p.pid() > 0) ? _h_pT_piplus : _h_pT_piminus)->fill(pT, 1.0/pT);
          ((p.pid() > 0) ? _tmp_pT_piplus : _tmp_pT_piminus)->fill(pT, 1.0/pT);
        }
      }

      const IdentifiedFinalState& protonfs = apply<IdentifiedFinalState>(event, "ProtonFS");
      for (const Particle& p : protonfs.particles()) {
        if (p.absrap() < 0.5) {
          /// @todo Use a binned counter to avoid this bin width cancellation hack
          const double pT = p.pT() / GeV;
          ((p.pid() > 0) ? _h_pT_proton : _h_pT_antiproton)->fill(pT, 1.0/pT);
        }
      }
      _sumWeightSelected->fill();
    }


    /// Finalize
    void finalize() {
      divide(_h_pT_piminus, _h_pT_piplus, _s_piminus_piplus);
      divide(_h_pT_antiproton, _h_pT_proton, _s_antipr_pr);
      divide(_h_pT_proton, _tmp_pT_piplus, _s_pr_piplus);
      divide(_h_pT_antiproton, _tmp_pT_piminus, _s_antipr_piminus);
      const YODA::Scatter1D factor = (1/(2*M_PI)) / *_sumWeightSelected;
      scale(_h_pT_piplus,     factor);
      scale(_h_pT_piminus,    factor);
      scale(_h_pT_proton,     factor);
      scale(_h_pT_antiproton, factor);
    }


  private:

    /// @{
    CounterPtr _sumWeightSelected;
    Histo1DPtr _h_pT_piplus, _h_pT_piminus, _h_pT_proton, _h_pT_antiproton;
    Histo1DPtr _tmp_pT_piplus, _tmp_pT_piminus;
    Scatter2DPtr _s_piminus_piplus, _s_antipr_pr, _s_pr_piplus, _s_antipr_piminus;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(STAR_2006_S6500200, STAR_2006_I709170);

}
