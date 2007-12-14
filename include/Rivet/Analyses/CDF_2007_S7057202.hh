// -*- C++ -*-
#ifndef RIVET_CDF_2007_S7057202_HH
#define RIVET_CDF_2007_S7057202_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"


namespace Rivet {

  /// @todo Needs full Doxygen commenting
  class CDF_2007_S7057202 : public Analysis {

  public:

    /// Constructor
    CDF_2007_S7057202():
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
      
      _ybins.push_back(0.1);
      _ybins.push_back(0.7);
      _ybins.push_back(1.1);
      _ybins.push_back(1.6);
      _ybins.push_back(2.1);
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
    
    /// Final state projection
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
        
    /// Rapidity bins of histograms
    /// @todo should be made constant (why isn't this being autobooked?)
    vector<double> _ybins;

    /// Counter for the number of events analysed (actually the sum of weights, hence double).
    double _eventsTried;
    
    /// The total generated cross section
    /// @todo Set the cross section from the generator
    double _xSecTot;
    
    // Histograms in different eta regions and the number of events
    // in each histogram
    //@{
    /// @todo Indexing a map by double is a bad idea
    map<double, AIDA::IHistogram1D*> _histosD07;
    map<AIDA::IHistogram1D*, double> _eventsPassedD07;

    AIDA::IHistogram1D* _histoD05;
    AIDA::IHistogram1D* _histoD10;

    double _eventsPassedD05;
    double _eventsPassedD10;
    //@}
  };

}

#endif
