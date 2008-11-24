// -*- C++ -*-
#ifndef RIVET_D0_2006_S6438750_HH
#define RIVET_D0_2006_S6438750_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/IsolationTools.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Inclusive isolated photon cross-section, differential in \f$ p_\perp(gamma) \f$.
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2006_S6438750 : public Analysis {

  public:

    /// Default constructor.
     D0_2006_S6438750();

    /// Factory method 
    static Analysis* create() {
      return new D0_2006_S6438750();
    }

    /// @name Publication metadata
    //@{
    /// Get a description of the analysis. 
    string spiresId() const {
      return "6438750";
    }
    /// Get a description of the analysis.
    string description() const {
      return "Inclusive isolated photon cross-section, differential in pT(gamma)";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2006";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("Phys.Lett.B639:151-158,2006, Erratum-ibid.B658:285-289,2008");
      ret.push_back("hep-ex/0511054 (plus erratum)");
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
    AIDA::IHistogram1D* _h_pTgamma;
    //@}

  };

}

#endif
