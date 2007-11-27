// -*- C++ -*-
#ifndef RIVET_ALEPH_1991_S2435284_HH
#define RIVET_ALEPH_1991_S2435284_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This analysis just measures the charged multiplicity.
  class ALEPH_1991_S2435284 : public Analysis {

  public:

    /// Default constructor.
    inline ALEPH_1991_S2435284()
      : _multproj(_fsproj)
    { 
      setBeams(ELECTRON, POSITRON); 
      addProjection(_fsproj);
      addProjection(_multproj);
    }

  public:


    /// Factory method
    static Analysis* create() { return new ALEPH_1991_S2435284(); }

    /// Return the name of this analysis.
    inline string getName() const {
      return "ALEPH_1991_S2435284";
    }

    virtual void init();

    virtual void analyze(const Event & event);

    virtual void finalize();

  private:

    /// @name The projectors used by this analysis.
    //{
    FinalState _fsproj;

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
