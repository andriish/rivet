// -*- C++ -*-
#ifndef RIVET_Run_HH
#define RIVET_Run_HH

#include "Rivet/AnalysisHandler.fhh"

namespace HepMC {
  class IO_GenEvent;
}

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
    Run& setListAnalyses(const bool& );
    
    bool prepareFile(const std::string& evtfile);
    bool processEvent(bool firstEvent);
    bool finalizeFile();
    
  private:

    /// AnalysisHandler object
    AnalysisHandler& _ah;
    
    double _xs;
    bool _listAnalyses;
    
    HepMC::IO_GenEvent* m_io;
    std::istream* m_istr;

  };


}

#endif
