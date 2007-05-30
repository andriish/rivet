// -*- C++ -*-
#ifndef RIVET_PL273B181_H
#define RIVET_PL273B181_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {

  /// This analysis just measures the charged multiplicity.
  class PL273B181 : public Analysis {

  public:

    /// Default constructor.
    inline PL273B181()
      : mult(fsproj)
    { 
      setBeams(ELECTRON, POSITRON); 
      addProjection(fsproj);
      addProjection(mult);
    }

  public:

    /// The name of this analysis is "PL273B181"
    inline string getName() const {
      return "PL273B181";
    }

    virtual void init();

    virtual void analyze(const Event & event);

    virtual void finalize();

  private:

    /// @name The projectors used by this analysis.
    //{
    FinalState fsproj;

    Multiplicity mult;
    //}

  private:

    /// Hide the assignment operator
    PL273B181 & operator=(const PL273B181 &);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* histChTot_;
    //@}

  };

}


#endif
