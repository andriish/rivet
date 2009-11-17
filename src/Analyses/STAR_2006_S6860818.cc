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


      IdentifiedFinalState k0sfs(-2.5, 2.5, 0.2*GeV);
      IdentifiedFinalState kaonfs(-2.5, 2.5, 0.2*GeV);
      IdentifiedFinalState lambdafs(-2.5, 2.5, 0.3*GeV);
      IdentifiedFinalState xifs(-2.5, 2.5, 0.5*GeV);
      IdentifiedFinalState omegafs(-2.5, 2.5, 0.5*GeV);
      k0sfs.acceptIdPair(K0S);
      kaonfs.acceptIdPair(KPLUS);
      lambdafs.acceptIdPair(LAMBDA);
      xifs.acceptIdPair(XIMINUS);
      omegafs.acceptIdPair(OMEGAMINUS);
      addProjection(k0sfs, "K0SFS");
      addProjection(kaonfs, "KAONFS");
      addProjection(lambdafs, "LAMBDAFS");
      addProjection(xifs, "XIFS");
      addProjection(omegafs, "OMEGAFS");

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

      const IdentifiedFinalState& k0sfs = applyProjection<IdentifiedFinalState>(event, "K0SFS");
      foreach (const Particle& p, k0sfs.particles()) {
        if (fabs(p.momentum().rapidity()) < 0.5) {
          const double pT = p.momentum().pT() / GeV;
          _h_pT_k0s->fill(pT, weight/pT);
        }
      }

      const IdentifiedFinalState& kaonfs = applyProjection<IdentifiedFinalState>(event, "KAONFS");
      foreach (const Particle& p, kaonfs.particles()) {
        if (fabs(p.momentum().rapidity()) < 0.5) {
          const double pT = p.momentum().pT() / GeV;
          if (p.pdgId()>0) {
            _h_pT_kplus->fill(pT, weight/pT);
          }
          else {
            _h_pT_kminus->fill(pT, weight/pT);
          }
        }
      }

      const IdentifiedFinalState& lambdafs = applyProjection<IdentifiedFinalState>(event, "LAMBDAFS");
      foreach (const Particle& p, lambdafs.particles()) {
        if (fabs(p.momentum().rapidity()) < 0.5) {
          const double pT = p.momentum().pT() / GeV;
          if (p.pdgId()>0) {
            _h_pT_lambda->fill(pT, weight/pT);
          }
          else {
            _h_pT_lambdabar->fill(pT, weight/pT);
          }
        }
      }

      const IdentifiedFinalState& xifs = applyProjection<IdentifiedFinalState>(event, "XIFS");
      foreach (const Particle& p, xifs.particles()) {
        if (fabs(p.momentum().rapidity()) < 0.5) {
          const double pT = p.momentum().pT() / GeV;
          if (p.pdgId()>0) {
            _h_pT_ximinus->fill(pT, weight/pT);
          }
          else {
            _h_pT_xiplus->fill(pT, weight/pT);
          }
        }
      }

      //const IdentifiedFinalState& omegafs = applyProjection<IdentifiedFinalState>(event, "OMEGAFS");
      //foreach (const Particle& p, omegafs.particles()) {
      //  if (fabs(p.momentum().rapidity()) < 0.5) {
      //    const double pT = p.momentum().pT() / GeV;
      //    _h_pT_omega->fill(pT, weight/pT);
      //  }
      //}

      _sumWeightSelected += event.weight();
    }


    /// Finalize
    void finalize() {
      AIDA::IHistogramFactory& hf = histogramFactory();
      const string dir = histoDir();

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
