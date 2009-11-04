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

    /// @name Set run properties
    //@{

    /// Get the cross-section for this run.
    Run& setCrossSection(const double xs);

    /// Declare whether to list available analyses
    Run& setListAnalyses(const bool dolist);

    //@}


    /// @name Get run conditions
    //@{

    /// Get beam IDs for this run, determined from first event
    const BeamPair& beams() const;

    /// Get energy for this run, determined from first event
    double sqrtS() const;

    //@}


    /// @name File processing stages
    //@{

    /// Set up HepMC file readers
    bool prepareFile(const std::string& evtfile);

    /// Handle next event
    bool processEvent(bool firstEvent);

    /// Close up HepMC I/O
    bool finalizeFile();

    //@}

    
  private:

    /// AnalysisHandler object
    AnalysisHandler& _ah;
    
    /// @name Run variables obtained from events or command line
    //@{

    /// Cross-section from command line
    double _xs;

    /// Centre of mass energy, determined from first event
    double _sqrts;

    /// Beam IDs, determined from first event
    BeamPair _beams;

    //@}


    /// Flag to show list of analyses
    bool _listAnalyses;


    /// @name HepMC I/O members
    //@{

    /// HepMC's own reader from streams   
    HepMC::IO_GenEvent* m_io;

    /// STL istream, used by IO_GenEvent if input is not a file
    std::istream* m_istr;

    //@}

  };


}

#endif
