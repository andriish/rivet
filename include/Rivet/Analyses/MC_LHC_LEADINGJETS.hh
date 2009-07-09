// -*- C++ -*-
#ifndef RIVET_MC_LHC_LEADINGJETS_HH
#define RIVET_MC_LHC_LEADINGJETS_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* Underlying event in leading jet, extended to the LHC
   * @author Andy Buckley
   */ 
  class MC_LHC_LEADINGJETS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_LHC_LEADINGJETS();

    /// Factory method
    static Analysis* create() {
      return new MC_LHC_LEADINGJETS();
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
    AIDA::IProfile1D *_hist_pcptsum;
    AIDA::IProfile1D *_hist_pmaxcptsum;
    AIDA::IProfile1D *_hist_pmincptsum;
    AIDA::IProfile1D *_hist_pcptave;

  };


}

#endif
