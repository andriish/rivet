// -*- C++ -*-
#ifndef RIVET_CDF_2007_S7057202_HH
#define RIVET_CDF_2007_S7057202_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/KtJets.hh"
#include "Rivet/Projections/FastJets.hh"


namespace Rivet {


  /// @todo Needs full Doxygen commenting
  class CDF_2007_S7057202 : public Analysis {

  public:

    /// Default constructor

    inline CDF_2007_S7057202():
    _fsproj(),
    _ktproj(_fsproj, KTTYPE, KTANGLE, KTRECOMBINATION, _ktRParam){
      setBeams(PROTON, ANTIPROTON);
      addProjection(_fsproj);
      addProjection(_ktproj);
      setNeedsCrossSection(true);
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
    
    FinalState _fsproj;
    KtJets _ktproj;
    
    //Parameters used in the KT algorithm
    enum KTParam {KTTYPE = 4, KTANGLE = 2, KTRECOMBINATION = 1};
    const static double _ktRParam;
    
    //Min jet PT cut
    const static double _jetMinPT;
    
    //Counter for the number of events analysed
    double _eventsTried;
    
    //The total generated cross section
    // @todo set the cross section from the generator
    double _xSecTot;
    
    //Histograms in different eta regions and the number of events
    //in each histogram
    
    map<double, AIDA::IHistogram1D*> _histos;
    map<AIDA::IHistogram1D*, double> _eventsPassed;

  };

}

#endif
