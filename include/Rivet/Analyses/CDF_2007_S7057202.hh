// -*- C++ -*-
#ifndef RIVET_CDF_2007_S7057202_HH
#define RIVET_CDF_2007_S7057202_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"


namespace Rivet {


  /// @todo Needs full Doxygen commenting
  class CDF_2007_S7057202 : public Analysis {

  public:

    /// Constructor.
    CDF_2007_S7057202():
      _ktRParam(0.7), _jetMinPT(54.0), 
      _fsproj(), _jetproj(_fsproj)
      /// @todo Reinstate KTTYPE, KTANGLE, KTRECOMBINATION, _ktRParam args to jet finder
      // This is a backwards use of an enum - completely unhelpful!
      //  enum KTParam { KTTYPE = 4, KTANGLE = 2, KTRECOMBINATION = 1 };
    {
      setBeams(PROTON, ANTIPROTON);
      addProjection(_fsproj);
      addProjection(_jetproj);
      setNeedsCrossSection(true);
    };
    
    /// Factory method.
    static Analysis* create() { 
      return new CDF_2007_S7057202(); 
    }

    /// Get the name of this analysis.
    string getName() const {
      return "CDF_2007_S7057202";
    }

  public:

    // Initializer.
    void init();

    // The analysis.
    void analyze(const Event& event);
    
    // Finalize histos.
    void finalize();

  private:

    /// Hide the assignment operator
    CDF_2007_S7057202& operator=(const CDF_2007_S7057202&);
    
    /// ...and the copy constructor
    CDF_2007_S7057202(const CDF_2007_S7057202&);

    // Parameters used in the KT algorithm.
    const double _ktRParam;
    
    /// Min jet \f$ p_T \f$ cut.
    const double _jetMinPT;
    
    /// Projections.
    FinalState _fsproj;
    FastJets _jetproj;
    
    /// Counter for the number of events analysed.
    double _eventsTried;
    
    /// The total generated cross section.
    /// @todo Set the cross section from the generator
    double _xSecTot;
    
    // Histograms in different eta regions and the number of events
    // in each histogram
    /// @todo Indexing a map by double is a bad idea...
    map<double, AIDA::IHistogram1D*> _histos;
    /// @todo Aaaargh!
    map<AIDA::IHistogram1D*, double> _eventsPassed;

  };

}

#endif
