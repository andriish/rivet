// -*- C++ -*-
#ifndef RIVET_PRD65092002_HH
#define RIVET_PRD65092002_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/TrackJet.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {

  class PRD65092002 : public Analysis {

  public:

    /// Default constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ and
    /// \f$ p_T > 0.5 \f$ GeV.
    inline PRD65092002()
      : _fsproj(-1.0, 1.0, 0.5), _trackjetproj(_fsproj) //, 
        //_dpsToward(0), _dpsAway(0), _dpsTrans(0), _numBins(50)
    { 
      setBeams(PROTON, ANTIPROTON);
      addProjection(_fsproj);
      addProjection(_trackjetproj);
    }


  public:

    /// Factory method
    static Analysis* create() { return new PRD65092002(); }

    /// Return the name of the analysis.
    inline string getName() const {
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

    /// @name Histogram collections
    //@{
    /// Profile histograms, binned in the \f$ p_T \f$ of the
    /// leading jet, for the pT sum in toward, transverse and away region
    AIDA::IProfile1D* _dataToward;
    AIDA::IProfile1D* _dataTrans;
    AIDA::IProfile1D* _dataAway;

    //@}


  private:

    /// Hide the assignment operator
    PRD65092002& operator=(const PRD65092002&);

  };

}

#endif
