// -*- C++ -*-
#ifndef RIVET_D0_2008_S7719523_HH
#define RIVET_D0_2008_S7719523_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/IsolationTools.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Measurement of isolated gamma + jet + X differential cross-sections
  /// Inclusive isolated gamma + jet cross-sections, differential in pT(gamma), for 
  /// various photon and jet rapidity bins.
  ///
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7719523 : public Analysis {

  public:

    /// Default constructor.
     D0_2008_S7719523();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7719523();
    }

    /// @name Publication metadata
    //@{
    /// Get a description of the analysis. 
    string getSpiresId() const {
      return "7719523";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Isolated gamma + jet cross-sections, differential in pT(gamma) for various y-bins";
    }    
    /// Experiment which performed and published this analysis. 
    string getExpt() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES). 
    string getYear() const {
      return "2008";
    }
    /// Publication references.
    vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("hep-ex/08081296");
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
    AIDA::IHistogram1D* _h_central_same_cross_section;
    AIDA::IHistogram1D* _h_central_opp_cross_section;
    AIDA::IHistogram1D* _h_forward_same_cross_section;
    AIDA::IHistogram1D* _h_forward_opp_cross_section;
    //@}

  };

}

#endif
