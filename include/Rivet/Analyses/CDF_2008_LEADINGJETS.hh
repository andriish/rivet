// -*- C++ -*-
#ifndef RIVET_CDF_2008_LEADINGJETS_HH
#define RIVET_CDF_2008_LEADINGJETS_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /* CDF Run II underlying event in leading jet events
   * @author Hendrik Hoeth
   * 
   * Rick Field's measurement of the underlying event in "leading jet" events.
   * The leading jet (CDF midpoint R=0.7) must be within |eta|<2 and defines
   * the "toward" phi direction. Particles are selected in |eta|<1. For the pT
   * related observables there is a pT>0.5 GeV cut. For sum(ET) there is no pT cut.
   * 
   * 
   * @par Run conditions
   * 
   * @arg \f$ \sqrt{s} = \f$ 1960 GeV
   * @arg Run with generic QCD events.
   * @arg Set particles with c*tau > 10 mm stable
   * @arg Several \f$ p_\perp^\text{min} \f$ cutoffs are probably required to fill the profile histograms:
   *   @arg \f$ p_\perp^\text{min} = \f$ 0 (min bias), 10, 20, 50, 100, 150 GeV
   *   @arg The corresponding merging points are at \f$ p_T = \f$ 0, 30, 50, 80, 130, 180 GeV
   * 
   */ 
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

    string name() const {
      return "CDF_2008_LEADINGJETS";
    }
    string description() const {
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

    AIDA::IProfile1D *_hist_pnchg;
    AIDA::IProfile1D *_hist_pmaxnchg;
    AIDA::IProfile1D *_hist_pminnchg;
    AIDA::IProfile1D *_hist_pdifnchg;
    AIDA::IProfile1D *_hist_pcptsum;
    AIDA::IProfile1D *_hist_pmaxcptsum;
    AIDA::IProfile1D *_hist_pmincptsum;
    AIDA::IProfile1D *_hist_pdifcptsum;
    AIDA::IProfile1D *_hist_pcptave;

  };


}

#endif
