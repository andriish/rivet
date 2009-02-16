// -*- C++ -*-
#ifndef RIVET_D0_2008_S7863608_HH
#define RIVET_D0_2008_S7863608_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/IsolationTools.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Measurement differential Z/gamma* + jet +X cross sections
  /// @author Gavin Hesketh
  class D0_2008_S7863608 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    D0_2008_S7863608();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7863608();
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "7863608";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Measurement of differential Z/gamma* + jet + X cross sections";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2008";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0808.1296 [hep-ex]");
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
    AIDA::IHistogram1D * _h_jet_pT_cross_section;
    AIDA::IHistogram1D * _h_jet_y_cross_section;
    AIDA::IHistogram1D * _h_Z_pT_cross_section;
    AIDA::IHistogram1D * _h_Z_y_cross_section;
    AIDA::IHistogram1D * _h_total_cross_section;
    //@}

  };


}

#endif
