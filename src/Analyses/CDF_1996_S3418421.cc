// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetYODA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF dijet angular distributions
  class CDF_1996_S3418421 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_1996_S3418421()
      : Analysis("CDF_1996_S3418421")
    {
    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(-4.2, 4.2);
      addProjection(FastJets(fs, FastJets::CDFJETCLU, 0.7), "Jets");

      _h_chi.addHistogram(241.0, 300.0, bookHisto1D(1, 1, 1));
      _h_chi.addHistogram(300.0, 400.0, bookHisto1D(1, 1, 2));
      _h_chi.addHistogram(400.0, 517.0, bookHisto1D(1, 1, 3));
      _h_chi.addHistogram(517.0, 625.0, bookHisto1D(1, 1, 4));
      _h_chi.addHistogram(625.0, 1800.0, bookHisto1D(1, 1, 5));

      _h_ratio = bookScatter2D(2, 1, 1);
      _chi_above_25.resize(_h_ratio->numPoints());
      _chi_below_25.resize(_h_ratio->numPoints());
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(50.0*GeV);
      if (jets.size()<2) {
        vetoEvent;
      }
      const FourMomentum jet1 = jets[0].momentum();
      const FourMomentum jet2 = jets[1].momentum();
      const double eta1 = jet1.eta();
      const double eta2 = jet2.eta();
      const double chi = exp(fabs(eta1 - eta2));
      if (fabs(eta2) > 2.0 || fabs(eta1) > 2.0 || chi > 5.0) vetoEvent;

      double m = FourMomentum(jet1 + jet2).mass();
      _h_chi.fill(m, chi, weight);

      // Fill ratio counter
      /// @todo Replace with a "binned counter": this is an inefficient bin lookup to avoid the "bin width factor"
      for (size_t i = 0; i < _h_ratio->numPoints(); ++i) {
        if (inRange(m/GeV, _h_ratio->point(i).xMin(), _h_ratio->point(i).xMax())) {
          ((chi > 2.5) ? _chi_above_25 : _chi_below_25)[i] += weight;
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      foreach (Histo1DPtr hist, _h_chi.getHistograms())
        normalize(hist);
      for (size_t bin = 0; bin < _h_ratio->numPoints(); ++bin)
        _h_ratio->point(bin).setY(_chi_below_25[bin]/_chi_above_25[bin]); //< @todo Add errors
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here
    std::vector<double> _chi_above_25;
    std::vector<double> _chi_below_25;

  private:

    /// @name Histograms
    //@{
    BinnedHistogram<double> _h_chi;
    Scatter2DPtr _h_ratio;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_1996_S3418421);

}
