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

    /// Default constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ and
    /// \f$ p_T > 0.5 \f$ GeV.
    inline PRD65092002()
      : _fsproj(-1.0, 1.0, 0.5), _trackjetproj(_fsproj), 
        _histToward(0), _histAway(0), _histTrans(0)
    { 
      setBeams(PROTON, ANTIPROTON);
    }


  public:

    /// The name of this analysis is "Test"
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

    //@{
    /// Histograms
    AIDA::IHistogram1D* _histToward;
    AIDA::IHistogram1D* _histAway;
    AIDA::IHistogram1D* _histTrans;
    //@}

    /// Internal storage of \f$ p_T \f$ flow data.
    class MiniHisto {
    public:
      /// Constructor
      inline MiniHisto() : numEntries(0), sumPt(0.0), sumPtSq(0.0) { }

      /// Appending operator
      inline void operator+=(const double evtPtSum) { 
        numEntries += 1;
        sumPt += evtPtSum;
        sumPtSq += evtPtSum*evtPtSum;
      }

      double numEntries;
      double sumPt;
      double sumPtSq;
    };

    /// These arrays of minimal histogram data, binned in the \f$ p_T \f$ of the
    /// leading jet, will be used to calculate the output profile histograms.
    //@{
    vector<MiniHisto> _dataToward;
    vector<MiniHisto> _dataAway;
    vector<MiniHisto> _dataTrans;
    //@}

  private:

    /// Hide the assignment operator
    PRD65092002& operator=(const PRD65092002& x);


  };

}

#endif
