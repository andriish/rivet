// -*- C++ -*-
#ifndef RIVET_ALEPH_1991_S2435284_HH
#define RIVET_ALEPH_1991_S2435284_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This analysis just measures the charged multiplicity.
  class ALEPH_1991_S2435284 : public Analysis {

  public:

    /// Default constructor.
    ALEPH_1991_S2435284()
      : _multproj(_cfsproj)
    { 
      setBeams(ELECTRON, POSITRON); 
      addProjection(_cfsproj);
      addProjection(_multproj);
    }

  public:


    /// Factory method
    static Analysis* create() { return new ALEPH_1991_S2435284(); }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "2435284";
    }
    /// Get a description of the analysis.
    //string getDescription() const {
    //  return "";
    //}
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "ALEPH";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
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

    /// @name The projectors used by this analysis.
    //{
    ChargedFinalState _cfsproj;
    Multiplicity _multproj;
    //}

  private:

    /// Hide the assignment operator
    ALEPH_1991_S2435284& operator=(const ALEPH_1991_S2435284&);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histChTot;
    //@}

  };

}


#endif
