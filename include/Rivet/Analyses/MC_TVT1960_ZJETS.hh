// -*- C++ -*-
#ifndef RIVET_MC_TVT1960_ZJETS_HH
#define RIVET_MC_TVT1960_ZJETS_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IsolationTools.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Monte Carlo validation observables for Z[e+ e-] + jets production at Tevatron Run II
  /// @author Frank Siegert
  class MC_TVT1960_ZJETS : public Analysis {

  public:

    /// Default constructor.
    MC_TVT1960_ZJETS();


    /// Factory method 
    static Analysis* create() {
      return new MC_TVT1960_ZJETS();
    }


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "NONE";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Monte Carlo validation observables for Z[e+ e-] + jets production at Tevatron Run II";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "MC";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "NONE";
    }
    /// Overwrite virtual name method for special analysis scheme
    string name() const {
      return "MC_TVT1960_ZJETS";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      return ret;
    }
    //@}


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
    AIDA::IHistogram1D * _h_jet1_pT;
    AIDA::IHistogram1D * _h_jet2_pT;
    AIDA::IHistogram1D * _h_jet3_pT;
    AIDA::IHistogram1D * _h_jet4_pT;
    AIDA::IHistogram1D * _h_jet_multi;
    AIDA::IHistogram1D * _h_deta_Z_jet1;
    AIDA::IHistogram1D * _h_dR_jet2_jet3;
    AIDA::IHistogram1D * _h_log10_d[4];
    AIDA::IDataPointSet * _h_log10_R[5];
    //@}

  };

}

#endif
