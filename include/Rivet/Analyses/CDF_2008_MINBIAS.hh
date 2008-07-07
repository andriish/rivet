// -*- C++ -*-
#ifndef RIVET_CDF_2008_MINBIAS_HH
#define RIVET_CDF_2008_MINBIAS_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief CDF Run II min-bias cross-section
  /// @author Hendrik Hoeth
  class CDF_2008_MINBIAS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.4 \f$ GeV.
    CDF_2008_MINBIAS()
    { 
      setBeams(PROTON, ANTIPROTON);
      const ChargedFinalState cfs(-1.0, 1.0, 0.4*GeV);
      addProjection(cfs, "FS");
    }


    /// Factory method
    static Analysis* create() {
      return new CDF_2008_MINBIAS();
    }
    //@}


  public:

    string getName() const {
      return "CDF_2008_MINBIAS";
    }
    string getDescription() const {
      return "CDF Run 2 min bias cross-section analysis";
    }

  public:

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    AIDA::IProfile1D *_hist_pt_vs_multiplicity;

  private:

    /// Hide the assignment operator.
    CDF_2008_MINBIAS& operator=(const CDF_2008_MINBIAS&);

  };

}

#endif
