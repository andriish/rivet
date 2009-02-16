// -*- C++ -*-
#ifndef RIVET_DELPHI_2002_069_CONF_603_HH
#define RIVET_DELPHI_2002_069_CONF_603_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief DELPHI b-fragmentation measurement
  /// @author Hendrik Hoeth
  class DELPHI_2002_069_CONF_603 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    DELPHI_2002_069_CONF_603() 
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "FS");
      addProjection(InitialQuarks(), "IQF");
    }


    /// Factory method.
    static Analysis* create() { 
      return new DELPHI_2002_069_CONF_603(); 
    }

    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string name() const {
      return "DELPHI_2002_069_CONF_603";
    }
    string spiresId() const {
      return "NONE";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Study of the b-quark fragmentation function at LEP I";
    }
    // Experiment which performed and published this analysis.
    string experiment() const {
     return "DELPHI";
    }
    // Collider on which the experiment was based.
    string collider() const {
     return "LEP1";
    }
    /// When published.
    string year() const {
     return "2002 (note)";
    }
    /// No journal or preprint references: this is a demo.
    vector<string> references() const {
      vector<string> ret;
      ret += "DELPHI note 2002-069-CONF-603 (ICHEP 2002)";
      return ret;
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the 
    /// inclusive single particle distributions' normalisations.

    AIDA::IHistogram1D *_histXbprim;
    AIDA::IHistogram1D *_histXbweak;

    AIDA::IProfile1D *_histMeanXbprim;
    AIDA::IProfile1D *_histMeanXbweak;

    //@}

  };

}

#endif
