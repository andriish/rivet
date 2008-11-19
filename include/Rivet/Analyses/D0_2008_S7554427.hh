// -*- C++ -*-
#ifndef RIVET_D0_2008_S7554427_HH
#define RIVET_D0_2008_S7554427_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"

namespace Rivet {


  /// @brief Measurement differntial Z/gamma* + jet +X cross sections
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7554427 : public Analysis {

  public:

    /// Default constructor.
    D0_2008_S7554427();


    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7554427();
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis. 
    string getSpiresId() const {
      return "7554427";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Measurement of differential Z/gamma* + X cross-section shape";
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
      ret.push_back("hep-ex/0712.0803");
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
    AIDA::IHistogram1D * _h_ZpT;
    AIDA::IHistogram1D * _h_forward_ZpT;
    //@}

  };


}

#endif
