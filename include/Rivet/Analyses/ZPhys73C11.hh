// -*- C++ -*-
#ifndef RIVET_ZPHYS73C11_HH
#define RIVET_ZPHYS73C11_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {


  class ZPhys73C11 : public Analysis {

  public:

    /// Default constructor.
    inline ZPhys73C11()
      : mult(fsproj), spher(fsproj) 
    { 
      setBeams(ELECTRON, POSITRON); 
      addProjection(fsproj);
      addProjection(mult);
      addProjection(spher);
    }

  public:

    /// Factory method
    static Analysis* create() { return new ZPhys73C11(); }

    /// The name of this analysis is "ZPhys73C11"
    inline string getName() const {
      return "ZPhys73C11";
    }

    virtual void init();

    virtual void analyze(const Event & event);

    virtual void finalize();

    /// Return the RivetInfo object of this analysis object.
    //virtual RivetInfo getInfo() const;

  private:

    /// The projectors used by this analysis.
    FinalState fsproj;

    Multiplicity mult;

    Sphericity spher;

  private:

    /// Hide the assignment operator
    ZPhys73C11 & operator=(const ZPhys73C11 &);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* histChTot_;
    AIDA::IHistogram1D* histSphericity_;
    AIDA::IHistogram1D* histPlanarity_;
    AIDA::IHistogram1D* histAplanarity_;
    //@}

  };

}


#endif
