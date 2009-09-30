// -*- C++ -*-
// CDF Run II inclusive jet cross-section using the midpoint algorithm.

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CDF Run II inclusive jet cross-section using the Midpoint algorithm.
  /// The analysis includes 1.1fb^-1 of CDF data and is the first with a 
  /// cone algorithm to include the forward region of the detector.
  /// arXiv:0807.2204 to be published in PRD
  class CDF_2008_S7828950 : public Analysis {
  public:
    
    /// Constructor
    CDF_2008_S7828950() : Analysis("CDF_2008_S7828950")
    {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
    }


    /// @name Analysis methods
    //@{

    // Book histos and set counters for number of events passed in each one
    void init() {
      const FinalState fs;
      addProjection(FastJets(fs, FastJets::CDFMIDPOINT, 0.7, 62.0*GeV), "JetsM07");

      /// @todo What actually are these histos showing?
      _binnedHistosR07.addHistogram(  0, 0.1, bookHistogram1D(1, 1, 1));
      _binnedHistosR07.addHistogram(0.1, 0.7, bookHistogram1D(2, 1, 1));
      _binnedHistosR07.addHistogram(0.7, 1.1, bookHistogram1D(3, 1, 1));
      _binnedHistosR07.addHistogram(1.1, 1.6, bookHistogram1D(4, 1, 1));
      _binnedHistosR07.addHistogram(1.6, 2.1, bookHistogram1D(5, 1, 1));

      size_t yind = 0;
      foreach (AIDA::IHistogram1D* hist, _binnedHistosR07.getHistograms()) {
        _yBinWidths[hist] = 2.0 * (_ybins[yind+1]-_ybins[yind]);
        ++yind;
      }
    }


    // Do the analysis
    void analyze(const Event& event) {
      const double weight = event.weight();    
      foreach (const Jet& jet, applyProjection<FastJets>(event, "JetsM07").jets()) {
        _binnedHistosR07.fill(fabs(jet.momentum().rapidity()), jet.momentum().pT(), weight);
      }
    }  


    // Normalise histograms to cross-section
    void finalize() {
      foreach (AIDA::IHistogram1D* hist, _binnedHistosR07.getHistograms()) {
        scale(hist, crossSection()/nanobarn/sumOfWeights()/_yBinWidths[hist]);
      }
    }

    //@}


  private:
    
    /// @name Histograms
    //@{

    /// The y bin width of each histogram
    map<AIDA::IHistogram1D*, double> _yBinWidths;

    /// The y bin edge values
    /// @todo Yuck!
    static const double _ybins[6];

    /// Histograms in different eta regions
    BinnedHistogram<double> _binnedHistosR07;

  };


  // Initialise static
  const double CDF_2008_S7828950::_ybins[] = { 0.0, 0.1, 0.7, 1.1, 1.6, 2.1 };

  
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_2008_S7828950> plugin_CDF_2008_S7828950;

}
