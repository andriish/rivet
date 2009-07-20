// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7828950_HH
#define RIVET_CDF_2008_S7828950_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// CDF Run II inclusive jet cross-section using the Midpoint algorithm.
  /// The analysis includes 1.1fb^-1 of CDF data and is the first with a 
  /// cone algorithm to include the forward region of the detector.
  /// arXiv:0807.2204 to be published in PRD
  class CDF_2008_S7828950 : public Analysis {
  public:
    
    /// @name Constructors etc.
    //@{
    
    /// Constructor
    CDF_2008_S7828950();
    
    /// Factory method
    static Analysis* create() { 
      return new CDF_2008_S7828950(); 
    }

    //@}
    

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}


  private:
    
    /// @name Histograms
    //@{

    /// The y bin width of each histogram
    map<AIDA::IHistogram1D*, double> _yBinWidths;

    /// The y bin edge values
    static const double _ybins[6];

    /// Histograms in different eta regions
    BinnedHistogram<double> _binnedHistosR07;

  };


}

#endif
