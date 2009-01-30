// -*- C++ -*-
#ifndef RIVET_ALEPH_1991_S2435284_HH
#define RIVET_ALEPH_1991_S2435284_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Multiplicity.hh"

namespace Rivet {


  /// @brief Measurement of ALEPH LEP1 charged multiplicity
  /// @author Andy Buckley
  class ALEPH_1991_S2435284 : public Analysis {
  public:

    /// Constructor.
    ALEPH_1991_S2435284() { 
      setBeams(ELECTRON, POSITRON); 
      const ChargedFinalState cfs;
      addProjection(cfs, "FS");
      addProjection(Multiplicity(cfs), "Mult");
    }

    /// Factory method.
    static Analysis* create() { 
      return new ALEPH_1991_S2435284(); 
    }

  public:

    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string spiresId() const {
      return "2435284";
    }
    /// Get a description of the analysis.
    string description() const {
      return "ALEPH LEP1 charged multiplicity measurement";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "ALEPH";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1991";
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event & event);
    virtual void finalize();
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histChTot;
    //@}

  };

}


#endif
