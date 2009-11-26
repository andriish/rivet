// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {

  /// @brief strange particle spectra in pp at 200 GeV
  class STAR_2006_S6860818 : public Analysis {
  public:

    /// Constructor
    STAR_2006_S6860818()
      : Analysis("STAR_2006_S6860818"),
        _sumWeightSelected(0.0)
    {
      setBeams(PROTON, PROTON);
    }

    /// Book projections and histograms
    void init() {
      ChargedFinalState bbc1(-5.0,-3.5, 0.0*GeV); // beam-beam-counter trigger
      ChargedFinalState bbc2( 3.5, 5.0, 0.0*GeV); // beam-beam-counter trigger
      addProjection(bbc1, "BBC1");
      addProjection(bbc2, "BBC2");

      UnstableFinalState ufs(-2.5, 2.5, 0.2*GeV);
      addProjection(ufs, "UFS");

      _h_pT_k0s        = bookHistogram1D(1, 1, 1);
      _h_pT_kminus     = bookHistogram1D(1, 2, 1);
      _h_pT_kplus      = bookHistogram1D(1, 3, 1);
      _h_pT_lambda     = bookHistogram1D(1, 4, 1);
      _h_pT_lambdabar  = bookHistogram1D(1, 5, 1);
      _h_pT_ximinus    = bookHistogram1D(1, 6, 1);
      _h_pT_xiplus     = bookHistogram1D(1, 7, 1);
      //_h_pT_omega      = bookHistogram1D(1, 8, 1);
    }


    /// Do the analysis
    void analyze(const Event& event) {
      const ChargedFinalState& bbc1 = applyProjection<ChargedFinalState>(event, "BBC1");
      const ChargedFinalState& bbc2 = applyProjection<ChargedFinalState>(event, "BBC2");
      if (bbc1.size()<1 || bbc2.size()<1) {
        getLog() << Log::DEBUG << "Failed beam-beam-counter trigger" << std::endl;
        vetoEvent;
      }

      const double weight = event.weight();

      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(event, "UFS");
      foreach (const Particle& p, ufs.particles()) {
        if (fabs(p.momentum().rapidity()) < 0.5) {
          const PdgId pid = p.pdgId();
          const double pT = p.momentum().pT() / GeV;
          switch (abs(pid)) {
            case K0S:
              if (pT > 0.2) {
                _h_pT_k0s->fill(pT, weight/pT);
              }
              break;
            case KPLUS:
              if (pT > 0.2) {
                pid > 0 ? _h_pT_kplus->fill(pT, weight/pT) : _h_pT_kminus->fill(pT, weight/pT);
              }
              break;
            case LAMBDA:
              if (pT > 0.3) {
                pid > 0 ? _h_pT_lambda->fill(pT, weight/pT) : _h_pT_lambdabar->fill(pT, weight/pT);
              }
              break;
            case XIMINUS:
              if (pT > 0.5) {
                pid > 0 ? _h_pT_ximinus->fill(pT, weight/pT) : _h_pT_xiplus->fill(pT, weight/pT);
              }
              break;
            //case OMEGAMINUS:
            //  if (pT > 0.5) {
            //    _h_pT_omega->fill(pT, weight/pT);
            //  }
            //  break;
          }

        }
      }

      _sumWeightSelected += event.weight();
    }


    /// Finalize
    void finalize() {
      //AIDA::IHistogramFactory& hf = histogramFactory();
      //const string dir = histoDir();
      //hf.divide(dir + "/d02-x01-y01", *_h_pT_piminus, *_h_pT_piplus);
      //hf.divide(dir + "/d02-x02-y01", *_h_pT_antiproton, *_h_pT_proton);
      //hf.divide(dir + "/d02-x03-y01", *_h_pT_proton, *_h_pT_piplus);
      //hf.divide(dir + "/d02-x04-y01", *_h_pT_antiproton, *_h_pT_piminus);

      scale(_h_pT_k0s,       1./(2*M_PI*_sumWeightSelected));
      scale(_h_pT_kminus,    1./(2*M_PI*_sumWeightSelected));
      scale(_h_pT_kplus,     1./(2*M_PI*_sumWeightSelected));
      scale(_h_pT_lambda,    1./(2*M_PI*_sumWeightSelected));
      scale(_h_pT_lambdabar, 1./(2*M_PI*_sumWeightSelected));
      scale(_h_pT_ximinus,   1./(2*M_PI*_sumWeightSelected));
      scale(_h_pT_xiplus,    1./(2*M_PI*_sumWeightSelected));
      //scale(_h_pT_omega,     1./(2*M_PI*_sumWeightSelected));
      getLog() << Log::DEBUG << "sumOfWeights()     = " << sumOfWeights() << std::endl;
      getLog() << Log::DEBUG << "_sumWeightSelected = " << _sumWeightSelected << std::endl;
    }

  private:

    double _sumWeightSelected;

    AIDA::IHistogram1D * _h_pT_k0s;
    AIDA::IHistogram1D * _h_pT_kminus;
    AIDA::IHistogram1D * _h_pT_kplus;
    AIDA::IHistogram1D * _h_pT_lambda;
    AIDA::IHistogram1D * _h_pT_lambdabar;
    AIDA::IHistogram1D * _h_pT_ximinus;
    AIDA::IHistogram1D * _h_pT_xiplus;
    //AIDA::IHistogram1D * _h_pT_omega;
  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<STAR_2006_S6860818> plugin_STAR_2006_S6860818;

}
