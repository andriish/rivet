// -*- C++ -*-
#ifndef RIVET_D0_2008_S7837160_HH
#define RIVET_D0_2008_S7837160_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"

namespace Rivet {


  /// @brief Measurement of W charge asymmetry from D0 Run II
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7837160 : public Analysis {

  public:

    /// Default constructor.
    D0_2008_S7837160();


    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7837160();
    }


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "7837160";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Measurement of W charge asymmetry from D0 Run II";
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
      //ret.push_back("Phys.Rev.Lett.");
      ret.push_back("arXiv:0807.3367");
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
    /// @todo Move into histo manager
    AIDA::IHistogram1D *_h_dsigplus_deta_25_35, *_h_dsigminus_deta_25_35;
    AIDA::IHistogram1D *_h_dsigplus_deta_35, *_h_dsigminus_deta_35;
    AIDA::IHistogram1D *_h_dsigplus_deta_25, *_h_dsigminus_deta_25;
    //@}

  };


}

#endif
