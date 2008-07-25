// -*- C++ -*-
#ifndef RIVET_CDF_2008_LEADINGJETS_HH
#define RIVET_CDF_2008_LEADINGJETS_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CDF Run II underlying event in leading jet events
  /// @author Hendrik Hoeth
  class CDF_2008_LEADINGJETS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    CDF_2008_LEADINGJETS()
    { 
      setBeams(PROTON, ANTIPROTON);

      // Final state for the jet finding
      const FinalState fsj(-4.0, 4.0, 0.0*GeV);
      addProjection(fsj, "FSJ");
      addProjection(FastJets(fsj, FastJets::CDFMIDPOINT, 0.7), "MidpointJets");

      // Final state for the sum(ET) distributions
      const FinalState fs(-1.0, 1.0, 0.0*GeV);
      addProjection(fs, "FS");

      // Charged final state for the distributions
      const ChargedFinalState cfs(-1.0, 1.0, 0.5*GeV);
      addProjection(cfs, "CFS");
    }


    /// Factory method
    static Analysis* create() {
      return new CDF_2008_LEADINGJETS();
    }
    //@}


  public:

    string getName() const {
      return "CDF_2008_LEADINGJETS";
    }
    string getDescription() const {
      return "CDF Run 2 underlying event in leading jet events";
    }

  public:

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    AIDA::IProfile1D *_hist_onchg;
    AIDA::IProfile1D *_hist_ocptsum;
    AIDA::IProfile1D *_hist_oetsum;

    //AIDA::IProfile1D *_hist_tnchg;
    //AIDA::IProfile1D *_hist_pnchg;
    //AIDA::IProfile1D *_hist_anchg;
    //AIDA::IProfile1D *_hist_tcptsum;
    //AIDA::IProfile1D *_hist_pcptsum;
    //AIDA::IProfile1D *_hist_acptsum;
    //AIDA::IProfile1D *_hist_tcptave;
    //AIDA::IProfile1D *_hist_pcptave;
    //AIDA::IProfile1D *_hist_acptave;
    //AIDA::IProfile1D *_hist_tcptmax;
    //AIDA::IProfile1D *_hist_pcptmax;
    //AIDA::IProfile1D *_hist_acptmax;

  private:

    /// Hide the assignment operator.
    CDF_2008_LEADINGJETS& operator=(const CDF_2008_LEADINGJETS&);

  };

}

#endif
