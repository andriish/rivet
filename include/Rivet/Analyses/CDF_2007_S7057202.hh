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
      _kttype(fastjet::kt_algorithm), _ktrecomb(fastjet::E_scheme),
      _ktRParam07(0.7), _ktRParam05(0.5), _ktRParam10(1.0), _jetMinPT(54.0),
      _fsproj(),
      _ktprojD07(_kttype, _ktrecomb, _ktRParam07, _fsproj),
      _ktprojD05(_kttype, _ktrecomb, _ktRParam05, _fsproj),
      _ktprojD10(_kttype, _ktrecomb, _ktRParam10, _fsproj) 
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
    
    /// Get the name of this analysis.
    inline string getName() const {
      return "CDF_2007_S7057202";
    }
    
  public:
    
    // Declaration of initialization method
    void init();
    
    // Declaration of analyzing method
    void analyze(const Event & event);
    
    void finalize();

  private:

    /// Hide the assignment operator
    CDF_2007_S7057202& operator=(const CDF_2007_S7057202&);
    
    /// ...and the copy constructor
    CDF_2007_S7057202(const CDF_2007_S7057202&);
    

    /// Parameters used in the KT algorithm('s)
    fastjet::JetFinder _kttype;
    fastjet::RecombinationScheme _ktrecomb;
    const double _ktRParam07;
    const double _ktRParam05;
    const double _ktRParam10;

    /// Min jet \f$ p_T \f$ cut.
    const double _jetMinPT;

    /// Projections.
    FinalState _fsproj;
    FastJets _ktprojD07;
    FastJets _ktprojD05;
    FastJets _ktprojD10;
        
    /// Rapidity bins of histograms
    /// @todo should be made constant
    vector<double> _ybins;


    /// Counter for the number of events analysed
    double _eventsTried;
    
    /// The total generated cross section
    /// @todo Set the cross section from the generator
    double _xSecTot;
    
    //Histograms in different eta regions and the number of events
    //in each histogram
    
    /// @todo Indexing a map by double is a bad idea
    map<double, AIDA::IHistogram1D*> _histosD07;
    /// @todo Aaaargh!
    map<AIDA::IHistogram1D*, double> _eventsPassedD07;

    AIDA::IHistogram1D* _histoD05;
    AIDA::IHistogram1D* _histoD10;

    double _eventsPassedD05;
    double _eventsPassedD10;

  };

}

#endif
