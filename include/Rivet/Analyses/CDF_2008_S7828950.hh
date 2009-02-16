// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7828950_HH
#define RIVET_CDF_2008_S7828950_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"


namespace Rivet {


  /// CDF Run II inclusive jet cross-section using the Midpoint algorithm.
  /// The analysis includes 1.1fb^-1 of CDF data and is the first with a 
  /// cone algorithm to include the forward region of the detector.
  /// arXiv:0807.2204 to be published in PRD
  /// @author Craig Group
  class CDF_2008_S7828950 : public Analysis {
  public:
    
    /// @name Constructors etc.
    //@{
    
    /// Constructor
    CDF_2008_S7828950()
      : _jetMinPT(62.0*GeV)
    {
      setBeams(PROTON, ANTIPROTON);
      //setSqrtS(1960*GeV);
      const FinalState fs;
      addProjection(FastJets(fs, FastJets::CDFMIDPOINT, 0.7), "JetsM07"); //??
      setNeedsCrossSection(true);
    }
    
    /// Factory method
    static Analysis* create() { 
      return new CDF_2008_S7828950(); 
    }
    //@}
    

    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "7828950";
    }
    /// A short description of the analysis.
    string summary() const {
      return "CDF Run II inclusive jet cross-section using the Midpoint algorithm";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2008";
    }
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0807.2204");
      //ret.push_back("PRD...");
      return ret;
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}


  private:

    /// Min jet \f$ p_T \f$ cut.
    /// @todo Make static const and UPPERCASE?
    const double _jetMinPT;
    
    /// Counter for the number of events analysed (actually the sum of weights, hence double).
    double _eventsTried;

    /// @name Histograms
    //@{
    /// The number of events in each histogram
    map<AIDA::IHistogram1D*, double> _eventsPassed;

    /// The y bin width of each histogram
    map<AIDA::IHistogram1D*, double> _yBinWidths;

    /// The y bin edge values
    static const double _ybins[6];

    /// Histograms in different eta regions
    BinnedHistogram<double> _binnedHistosR07;

  };


}

#endif
