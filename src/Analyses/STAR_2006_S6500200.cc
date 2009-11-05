// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {

  /// @brief identified hadron spectra in pp at 200 GeV
  class STAR_2006_S6500200 : public Analysis {
  public:

    /// Constructor
    STAR_2006_S6500200()
      : Analysis("STAR_2006_S6500200")
    {
      setBeams(PROTON, PROTON);
    }

    /// Book projections and histograms
    void init() {
      IdentifiedFinalState pionfs(-3.0, 3.0, 0.3*GeV);
      IdentifiedFinalState protonfs(-3.0, 3.0, 0.4*GeV);
      pionfs.acceptIdPair(PIPLUS);
      protonfs.acceptIdPair(PROTON);
      addProjection(pionfs, "PIONFS");
      addProjection(protonfs, "PROTONFS");

      _h_pT_piplus     = bookHistogram1D(1, 1, 1);
      _h_pT_piminus    = bookHistogram1D(1, 2, 1);
      _h_pT_proton     = bookHistogram1D(1, 3, 1);
      _h_pT_antiproton = bookHistogram1D(1, 4, 1);
    }


    /// Do the analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const IdentifiedFinalState& pionfs = applyProjection<IdentifiedFinalState>(event, "PIONFS");
      const IdentifiedFinalState& protonfs = applyProjection<IdentifiedFinalState>(event, "PROTONFS");
      foreach (const Particle& p, pionfs.particles()) {
        if (fabs(p.momentum().rapidity()) < 0.5) {
          const double pT = p.momentum().pT() / GeV;
          if (p.pdgId()>0) {
            _h_pT_piplus->fill(pT, weight/pT);
          }
          else {
            _h_pT_piminus->fill(pT, weight/pT);
          }
        }
      }
      foreach (const Particle& p, protonfs.particles()) {
        if (fabs(p.momentum().rapidity()) < 0.5) {
          const double pT = p.momentum().pT() / GeV;
          if (p.pdgId()>0) {
            _h_pT_proton->fill(pT, weight/pT);
          }
          else {
            _h_pT_antiproton->fill(pT, weight/pT);
          }
        }
      }
    }


    /// Finalize
    void finalize() {
      AIDA::IHistogramFactory& hf = histogramFactory();
      const string dir = histoDir();

      hf.divide(dir + "/d02-x01-y01", *_h_pT_piminus, *_h_pT_piplus);
      hf.divide(dir + "/d02-x02-y01", *_h_pT_antiproton, *_h_pT_proton);
      hf.divide(dir + "/d02-x03-y01", *_h_pT_proton, *_h_pT_piplus);
      hf.divide(dir + "/d02-x04-y01", *_h_pT_antiproton, *_h_pT_piminus);

      scale(_h_pT_piplus,     1./(2*M_PI*sumOfWeights()));
      scale(_h_pT_piminus,    1./(2*M_PI*sumOfWeights()));
      scale(_h_pT_proton,     1./(2*M_PI*sumOfWeights()));
      scale(_h_pT_antiproton, 1./(2*M_PI*sumOfWeights()));
    }

  private:

    AIDA::IHistogram1D * _h_pT_piplus;
    AIDA::IHistogram1D * _h_pT_piminus;
    AIDA::IHistogram1D * _h_pT_proton;
    AIDA::IHistogram1D * _h_pT_antiproton;
  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<STAR_2006_S6500200> plugin_STAR_2006_S6500200;

}
