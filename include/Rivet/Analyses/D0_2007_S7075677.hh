// -*- C++ -*-
#ifndef RIVET_D0_2007_S7075677_HH
#define RIVET_D0_2007_S7075677_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"

namespace Rivet {


  /// @brief Measurement of D0 Run II Z pT diff cross-section shape
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2007_S7075677 : public Analysis {

  public:

    /// Default constructor.
    D0_2007_S7075677();


    /// Factory method 
    static Analysis* create() {
      return new D0_2007_S7075677();
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis. 
    string spiresId() const {
      return "7075677";
    }
    /// Get a description of the analysis.
    string description() const {
      return "Z/gamma* + X cross-section shape, differential in y(Z)";
    }
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2007";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("Phys.Rev.D76:012003,2007");
      ret.push_back("hep-ex/0702025");
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
    AIDA::IHistogram1D * _h_yZ;
    //@}

  };


}

#endif
