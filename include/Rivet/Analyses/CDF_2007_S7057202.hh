// -*- C++ -*-
#ifndef RIVET_CDF_2007_S7057202_HH
#define RIVET_CDF_2007_S7057202_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// CDF Run II inclusive jet cross-section using the kT algorithm.
  /// @author James Monk
  class CDF_2007_S7057202 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2007_S7057202()
      : _minY(0.1), _maxY(0.7), _jetMinPT(54.0*GeV)
   {
      setBeams(PROTON, ANTIPROTON);
      //setSqrtS(1960*GeV);
      const FinalState fs;
      addProjection(FastJets(fs, FastJets::KT, 0.5), "JetsD05");
      addProjection(FastJets(fs, FastJets::KT, 0.7), "JetsD07");
      addProjection(FastJets(fs, FastJets::KT, 1.0), "JetsD10");
      setNeedsCrossSection(true);
    }
    
    /// Factory method
    static Analysis* create() { 
      return new CDF_2007_S7057202(); 
    }
    //@}
    

    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "7057202";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "CDF Run II inclusive jet cross-section using the kT algorithm";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2007";
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}


  private:

    /// Rapidity range of histograms for R=0.05 and R=1 kt jets
    const double _minY, _maxY;
        
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
    BinnedHistogram<double> _binnedHistosD07;

    // Single histogram for the \f$R=0.5\f$ \f$k_\perp\f$ jets
    AIDA::IHistogram1D* _histoD05;

    // Single histogram for the \f$R=1.0\f$ \f$k_\perp\f$ jets
    AIDA::IHistogram1D* _histoD10;
    //@}

  private:
    /// Hide the assignment operator
    CDF_2007_S7057202& operator=(const CDF_2007_S7057202&);
  };

}

#endif
