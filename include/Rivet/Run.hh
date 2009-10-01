// -*- C++ -*-
#ifndef RIVET_Run_HH
#define RIVET_Run_HH

#include "Rivet/AnalysisHandler.fhh"

namespace Rivet {


  class Run {

  public:

    /// @name Standard constructors and destructors. */
    //@{
    /// The standard constructor.
    Run(AnalysisHandler& ah);

    /// The destructor
    ~Run();
    //@}


  public:

    /// Get the name of this run.
    Run& setCrossSection(const double& xs);
    Run& setMaxEvtNum(const int& n);
    Run& setListAnalyses(const bool& );
    
    bool processFile(const std::string& evtfile);
    
  private:
    
    void logNEvt();


  private:

    /// AnalysisHandler object
    AnalysisHandler& _ah;
    
    double _xs;
    long _maxEvtNum, _numEvents;
    bool _listAnalyses;

  };


}

#endif
