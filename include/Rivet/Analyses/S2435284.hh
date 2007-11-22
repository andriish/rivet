// -*- C++ -*-
#ifndef RIVET_S2435284_HH
#define RIVET_S2435284_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This analysis just measures the charged multiplicity.
  class S2435284 : public Analysis {

  public:

    /// Default constructor.
    inline S2435284()
      : _multproj(_fsproj)
    { 
      setBeams(ELECTRON, POSITRON); 
      addProjection(_fsproj);
      addProjection(_multproj);
    }

  public:


    /// Factory method
    static Analysis* create() { return new S2435284(); }

    /// Return the name of this analysis.
    inline string getName() const {
      return "S2435284";
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
    S2435284& operator=(const S2435284&);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histChTot;
    //@}

  };

}


#endif
