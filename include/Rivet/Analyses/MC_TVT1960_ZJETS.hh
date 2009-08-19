// -*- C++ -*-
#ifndef RIVET_MC_TVT1960_ZJETS_HH
#define RIVET_MC_TVT1960_ZJETS_HH

#include "Rivet/Analyses/MC_JetAnalysis.hh"

namespace Rivet {


  /// @brief Monte Carlo validation observables for Z[e+ e-] + jets production at Tevatron Run II
  /// @author Frank Siegert
  class MC_TVT1960_ZJETS : public MC_JetAnalysis {

  public:

    /// Default constructor.
    MC_TVT1960_ZJETS();


    /// Factory method 
    static Analysis* create() {
      return new MC_TVT1960_ZJETS();
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
    AIDA::IHistogram1D * _h_Z_mass;
    AIDA::IHistogram1D * _h_Z_pT;
    AIDA::IHistogram1D * _h_Z_y;
    AIDA::IHistogram1D * _h_Z_jet1_deta;
    AIDA::IHistogram1D * _h_Z_jet1_dR;
    AIDA::IHistogram1D * _h_lepton_pT;
    AIDA::IHistogram1D * _h_lepton_eta;
    //@}

  };

}

#endif
