// -*- C++ -*-
#ifndef RIVET_HepEx0701051_HH
#define RIVET_HepEx0701051_HH

#include "Rivet/Analysis.hh"

#include "Rivet/Projections/KtJets.hh"
#ifdef HAVE_FASTJET	
#include "Rivet/Projections/FastJets.hh"
#endif


namespace Rivet {


  /// @todo Needs full Doxygen commenting
  class HepEx0701051 : public Analysis {

  public:

    /// Default constructor

    inline HepEx0701051():
    _fsproj(),
    _ktproj(_fsproj, KTTYPE, KTANGLE, KTRECOMBINATION, _ktRParam){
      setBeams(PROTON, ANTIPROTON);
      addProjection(_fsproj);
      addProjection(_ktproj);
    };
    
    /// Factory method
    static Analysis* create() { 
      return new HepEx0701051(); 
    }

    /// Get the name of this analysis.
    inline string getName() const {
      return "HepEx0701051";
    }

  public:

    // Declaration of initialization method
    void init();

    // Declaration of analyzing method
    void analyze(const Event & event);
    
    void finalize();

  private:

    /// Hide the assignment operator
    HepEx0701051& operator=(const HepEx0701051&);
    
    /// ...and the copy constructor
    HepEx0701051(const HepEx0701051&);
    
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
