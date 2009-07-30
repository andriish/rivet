// -*- C++ -*-
#ifndef RIVET_MC_LHC_ZANALYSIS_HH
#define RIVET_MC_LHC_ZANALYSIS_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    
  class MC_LHC_ZANALYSIS : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// constructor
    MC_LHC_ZANALYSIS();

    /// Factory method
    static Analysis* create() { 
      return new MC_LHC_ZANALYSIS(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

private:
AIDA::IHistogram1D* _hist_chargemulti;
AIDA::IHistogram1D* _hist_chargept;
AIDA::IHistogram1D* _hist_chargemeanpt;
AIDA::IHistogram1D* _hist_chargermspt;
AIDA::IHistogram1D* _hist_zcount;
AIDA::IHistogram1D* _hist_zpt;
AIDA::IHistogram1D* _hist_zlogpt;
//AIDA::IHistogram1D* _hist_zpthigh;
// AIDA::IHistogram1D* _hist_zlogpthigh;
AIDA::IHistogram1D* _hist_zeta;
AIDA::IHistogram1D* _hist_zphi;
AIDA::IHistogram1D* _hist_zmass;
AIDA::IHistogram1D* _hist_zlogmass;
AIDA::IHistogram1D* _hist_jetcount;
AIDA::IHistogram1D* _hist_jetpt;
AIDA::IHistogram1D* _hist_jetlogpt;

  };

}

#endif

