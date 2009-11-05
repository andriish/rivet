// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CDF_1996_S3418421 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_1996_S3418421()
      : Analysis("CDF_1996_S3418421"), _nevt(0)
    {
      setBeams(PROTON, ANTIPROTON);
    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(-4.2, 4.2);
      addProjection(FastJets(fs, FastJets::CDFJETCLU, 0.7, 0.0*GeV), "Jets");

      _h_chi.addHistogram(241.0, 300.0, bookHistogram1D(1, 1, 1));
      _h_chi.addHistogram(300.0, 400.0, bookHistogram1D(1, 1, 2));
      _h_chi.addHistogram(400.0, 517.0, bookHistogram1D(1, 1, 3));
      _h_chi.addHistogram(517.0, 625.0, bookHistogram1D(1, 1, 4));
      _h_chi.addHistogram(625.0, 1800.0, bookHistogram1D(1, 1, 5));
      
      _h_ratio = bookHistogram1D(2,1,1);
      _chi_above_25.resize(_h_ratio->axis().bins());
      _chi_below_25.resize(_h_ratio->axis().bins());
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      Jets jets = applyProjection<FastJets>(event, "Jets").jetsByEt(50.0*GeV);
      if (jets.size()<2) {
        vetoEvent;
      }
      FourMomentum jet1 = jets[0].momentum();
      FourMomentum jet2 = jets[1].momentum();
      double eta1 = fabs(jet1.eta());
      double eta2 = fabs(jet2.eta());
      double chi = exp(fabs(eta1-eta2));
      if (eta2>2.0 || eta1>2.0 || chi>5.0) {
        vetoEvent;
      }
      
      double m = FourMomentum(jet1+jet2).mass();
      _h_chi.fill(m, chi, weight);
      
      // fill ratio counter
      int bin = _h_ratio->coordToIndex(m);
      _nevt++;
      if (chi>2.5) {
        _chi_above_25[bin] += weight;
      }
      else {
        _chi_below_25[bin] += weight;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      
      foreach (AIDA::IHistogram1D* hist, _h_chi.getHistograms()) {
        // because HepData contains 100/N instead of 1/N this is 100 in the aida
        normalize(hist, 100.0);
      }

      // now mimic filling of the ratio histogram
      std::vector<double> ratios(_h_ratio->axis().bins());
      for (size_t i=0; i<ratios.size(); ++i) {
        ratios[i] = _chi_below_25[i]/_chi_above_25[i]/double(_nevt);
      }
      for (size_t bin=0; bin<ratios.size(); ++bin) {
        double coord = _h_ratio->binMean(bin);
        for (size_t n=0; n<_nevt; ++n) {
          _h_ratio->fill(coord, ratios[bin]);
        }
      }
      
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here
    std::vector<double> _chi_above_25;
    std::vector<double> _chi_below_25;
    size_t _nevt;

  private:

    /// @name Histograms
    //@{
    BinnedHistogram<double> _h_chi;
    AIDA::IHistogram1D* _h_ratio;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_1996_S3418421> plugin_CDF_1996_S3418421;


}
