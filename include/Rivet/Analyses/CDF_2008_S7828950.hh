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

 /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Measurement of the inclusive jet cross section in $p\bar{p}$ "
         << "collisions at $\sqrt{s}=1.96$ TeV as a function of jet Et, "
         <<  "for Et $>$ 62 GeV. "
         << "The data is collected by "
	 << "the CDF II detector and has an integrated luminosity of "
         << "1.13 fb$^{-1}$. The measurement was made using the cone-based "
         << "Midpoint jet clustering algorithm in rapidity bins within "
         <<  "$|y|<2.1$. This measurement can be used to provide increased "
         << "precision in PDFs at high parton momentum fraction $x$."; 
      return os.str();
    }
  /// Characteristics of events to be processed by this analysis
    string runInfo() const {
      ostringstream os;
      os << "Requires 2->2 QCD scattering processes. "
         << "The minimum jet Et is 62 GeV, so a cut on kinematic pTmin "
         << "may be required for good statistics.";
      return os.str();
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 2";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2008";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Craig Group <group@fnal.gov>";
      return rtn;
    }

    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0807.2204");
      ret.push_back("Phys.Rev.D78:052006,2008");
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
