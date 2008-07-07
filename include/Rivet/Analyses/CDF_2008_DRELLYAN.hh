// -*- C++ -*-
#ifndef RIVET_CDF_2008_DRELLYAN_HH
#define RIVET_CDF_2008_DRELLYAN_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {


  /// @brief CDF Run II underlying event in Drell-Yan
  /// @author Hendrik Hoeth
  class CDF_2008_DRELLYAN : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.5 \f$ GeV.
    CDF_2008_DRELLYAN()
    { 
      setBeams(PROTON, ANTIPROTON);
      const ChargedFinalState cfs(-1.0, 1.0, 0.5*GeV);
      addProjection(cfs, "FS");
      addProjection(ChargedLeptons(cfs), "CL");
    }


    /// Factory method
    static Analysis* create() {
      return new CDF_2008_DRELLYAN();
    }
    //@}


  public:

    string getName() const {
      return "CDF_2008_DRELLYAN";
    }
    string getDescription() const {
      return "CDF Run 2 underlying event in Drell-Yan";
    }

  public:

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    AIDA::IProfile1D *_hist_tnchg;
    AIDA::IProfile1D *_hist_pnchg;
    AIDA::IProfile1D *_hist_anchg;
    AIDA::IProfile1D *_hist_tcptsum;
    AIDA::IProfile1D *_hist_pcptsum;
    AIDA::IProfile1D *_hist_acptsum;
    AIDA::IProfile1D *_hist_tcptave;
    AIDA::IProfile1D *_hist_pcptave;
    AIDA::IProfile1D *_hist_acptave;
    AIDA::IProfile1D *_hist_tcptmax;
    AIDA::IProfile1D *_hist_pcptmax;
    AIDA::IProfile1D *_hist_acptmax;

  private:

    /// Hide the assignment operator.
    CDF_2008_DRELLYAN& operator=(const CDF_2008_DRELLYAN&);

  };

}

#endif
