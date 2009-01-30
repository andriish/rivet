// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7540469_HH
#define RIVET_CDF_2008_S7540469_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/IsolationTools.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Measurement differential Z/gamma* + jet + X cross sections
  /// @author Frank Siegert
  class CDF_2008_S7540469 : public Analysis {

  public:

    /// Default constructor.
    CDF_2008_S7540469();


    /// Factory method 
    static Analysis* create() {
      return new CDF_2008_S7540469();
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis. 
    string spiresId() const {
      return "7540469";
    }
    /// Get a description of the analysis.
    string description() const {
      return "Measurement of differential Z/gamma* + jet + X cross sections";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2008";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0711.3717 [hep-ex]");
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
    AIDA::IHistogram1D * _h_jet_multiplicity;
    AIDA::IHistogram1D * _h_jet_pT_cross_section_incl_1jet;
    AIDA::IHistogram1D * _h_jet_pT_cross_section_incl_2jet;
    //@}

  };

}

#endif
