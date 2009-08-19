// -*- C++ -*-
#ifndef RIVET_MC_TVT1960_PHOTONJETS_HH
#define RIVET_MC_TVT1960_PHOTONJETS_HH

#include "Rivet/Analyses/MC_JetAnalysis.hh"

namespace Rivet {


  class MC_TVT1960_PHOTONJETS : public MC_JetAnalysis {

  public:

    /// Default constructor.
    MC_TVT1960_PHOTONJETS();


    /// Factory method 
    static Analysis* create() {
      return new MC_TVT1960_PHOTONJETS();
    }


    /// @name Analysis methods
    //@{ 
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_photon_pT;
    AIDA::IHistogram1D * _h_photon_y;
    AIDA::IHistogram1D * _h_photon_jet1_deta;
    AIDA::IHistogram1D * _h_photon_jet1_dR;
    //@}

  };

}

#endif
