// -*- C++ -*-
#ifndef RIVET_PL273B181_H
#define RIVET_PL273B181_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/RivetAIDA.hh"
#include "AIDA/IHistogram1D.h"



namespace Rivet {

  /// This class just measures the charged multiplicity
  class PL273B181 : public Analysis {

  public:

    /// Default constructor.
    inline PL273B181();

    /// Copy constructor.
    inline PL273B181(const PL273B181 &);

    /// Destructor
    virtual ~PL273B181();

  public:

    /// The name of this analysis is "PL273B181"
    inline std::string name() const {
      return "PL273B181";
    }

    virtual void init();

    virtual void analyze(const Event & event);

    virtual void finalize();

    /// Return the RivetInfo object of this analysis object.
    virtual RivetInfo getInfo() const;

  private:

    /// The projectors used by this analysis.
    Multiplicity mult;

    Sphericity spher;

  private:

    // Hide the assignment operator
    PL273B181 & operator=(const PL273B181 &);

    // Histograms
    AIDA::IHistogram1D* histChTot_;
    AIDA::IHistogram1D* histSphericity_;
    AIDA::IHistogram1D* histPlanarity_;
    AIDA::IHistogram1D* histAplanarity_;

  };

}

#include "Rivet/Analysis/PL273B181.icc"

#endif
