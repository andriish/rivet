// -*- C++ -*-
#ifndef RIVET_CDF_2007_S7057202_HH
#define RIVET_CDF_2007_S7057202_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  /// @todo Needs full Doxygen commenting
  class CDF_2007_S7057202 : public Analysis {

  public:

    /// Constructor
    CDF_2007_S7057202():
      _minY(0.1), _maxY(0.7),
      _fsproj(),
      _ktprojD07(_fsproj, FastJets::KT, 0.7),
      _ktprojD05(_fsproj, FastJets::KT, 0.5),
      _ktprojD10(_fsproj, FastJets::KT, 1.0),
      _jetMinPT(54.0*GeV)
   {
      setBeams(PROTON, ANTIPROTON);
      addProjection(_fsproj);
      addProjection(_ktprojD07);
      addProjection(_ktprojD05);
      addProjection(_ktprojD10);
      setNeedsCrossSection(true);
    };

    
    /// Factory method
    static Analysis* create() { 
      return new CDF_2007_S7057202(); 
    }
    

    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "7057202";
    }
    /// Get a description of the analysis.
    //string getDescription() const {
    //  return "";
    //}
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

    /// Hide the assignment operator
    CDF_2007_S7057202& operator=(const CDF_2007_S7057202&);

    /// Rapidity range of histograms for R=0.05 and R=1 kt jets
    const double _minY, _maxY;
    
    /// Projections.
    FinalState _fsproj;
    /// @name Jet projections (with different R)
    //@{
    FastJets _ktprojD07;
    FastJets _ktprojD05;
    FastJets _ktprojD10;
    //@}
    
    /// Min jet \f$ p_T \f$ cut.
    /// @todo Make static const and UPPERCASE?
    const double _jetMinPT;
    
    /// Counter for the number of events analysed (actually the sum of weights, hence double).
    double _eventsTried;
    ///The number of events in each histogram
    map<AIDA::IHistogram1D*, double> _eventsPassed;
    ///Histograms in different eta regions
    BinnedHistogram<double> _binnedHistosD07;
    //@{Single histograms for the R=0.5 and R=1.0 KT jets
    AIDA::IHistogram1D* _histoD05;
    AIDA::IHistogram1D* _histoD10;
    //@}
  };

}

#endif
