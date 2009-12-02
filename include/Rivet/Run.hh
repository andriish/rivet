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
    bool init(const std::string& evtfile);

    /// Read the next HepMC event
    bool readEvent();

    /// Handle next event
    bool processEvent();

    /// Close up HepMC I/O
    bool finalize();

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

    /// Current event
    shared_ptr<GenEvent> _evt;

    /// Output stream for HepMC writer
    shared_ptr<std::istream> _istr;

    /// HepMC I/O writer
    shared_ptr<HepMC::IO_GenEvent> _io;

    //@}

  };


}

#endif
