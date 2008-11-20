// -*- C++ -*-
#ifndef RIVET_D0_2008_S7662670_HH
#define RIVET_D0_2008_S7662670_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/IsolationTools.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Measurement of D0 differential jet cross sections
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7662670 : public Analysis {

  public:

    /// Default constructor.
    D0_2008_S7662670();


    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7662670();
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis. 
    string getSpiresId() const {
      return "7662670";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Measurement of D0 Run II differential jet cross sections";
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
      ret.push_back("arXiv:0802.2400");
      ret.push_back("Phys.Rev.Lett.101:062001,2008");
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
    AIDA::IHistogram1D * _h_dsigdptdy_y00_04;
    AIDA::IHistogram1D * _h_dsigdptdy_y04_08;
    AIDA::IHistogram1D * _h_dsigdptdy_y08_12;
    AIDA::IHistogram1D * _h_dsigdptdy_y12_16;
    AIDA::IHistogram1D * _h_dsigdptdy_y16_20;
    AIDA::IHistogram1D * _h_dsigdptdy_y20_24;
    //@}

  };

}

#endif
