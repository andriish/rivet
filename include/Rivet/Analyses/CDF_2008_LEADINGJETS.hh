// -*- C++ -*-
#ifndef RIVET_CDF_2008_LEADINGJETS_HH
#define RIVET_CDF_2008_LEADINGJETS_HH

#include "Rivet/Analysis.hh"

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

    CDF_2008_LEADINGJETS();

    /// Factory method
    static Analysis* create() {
      return new CDF_2008_LEADINGJETS();
    }
    //@}


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
