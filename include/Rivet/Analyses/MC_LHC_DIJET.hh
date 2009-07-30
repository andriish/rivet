// -*- C++ -*-
#ifndef RIVET_MC_LHC_DIJET_HH
#define RIVET_MC_LHC_DIJET_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    
  class MC_LHC_DIJET : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// constructor
    MC_LHC_DIJET();

    /// Factory method
    static Analysis* create() { 
      return new MC_LHC_DIJET(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

private:
AIDA::IHistogram1D* _hist_jetcount;
AIDA::IHistogram1D* _hist_jetpt;
AIDA::IHistogram1D* _hist_jetptlog;
AIDA::IHistogram1D* _hist_leadingjetpt;
AIDA::IHistogram1D* _hist_secondleadingjetpt;
AIDA::IHistogram1D* _hist_jetphi;
AIDA::IHistogram1D* _hist_jetdphi;
AIDA::IHistogram1D* _hist_jeteta;
AIDA::IHistogram1D* _hist_jetdeta;
AIDA::IHistogram1D* _hist_chargemultiplicity;
AIDA::IHistogram1D* _hist_chargemeanpt;
AIDA::IHistogram1D* _hist_chargept;
AIDA::IHistogram1D* _hist_chargelogpt;
AIDA::IHistogram1D* _hist_chargermspt;

  };

}

#endif

