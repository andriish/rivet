// -*- C++ -*-
#ifndef RIVET_PRD65092002_HH
#define RIVET_PRD65092002_HH

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/TrackJet.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {

  class PRD65092002 : public Analysis {

  public:

    /// Default constructor
    inline PRD65092002()
      : _fsproj(-1.0, 1.0, 0.5), _trackjetproj(_fsproj), 
        _histToward(0), _histAway(0), _histTrans(0)
    { }

    // cuts as args
    // fsproj: -1 < eta < 1 and PT > 0.5 GeV

  public:

    /// The name of this analysis is "Test"
    inline string name() const {
      return "PRD65092002";
    }

  public:

    void init();
    
    void analyze(const Event & event);
    
    void finalize();


  private:

    /// The FinalState projection used.
    FinalState _fsproj;

    /// The TrackJet projection used by this analysis.
    TrackJet _trackjetproj;


  private:

    //@{
    /// Histograms
    AIDA::IHistogram1D* _histToward;
    AIDA::IHistogram1D* _histAway;
    AIDA::IHistogram1D* _histTrans;
    //@}


  private:

    /// Hide the assignment operator
    PRD65092002& operator=(const PRD65092002& x);


  };

}

#endif
