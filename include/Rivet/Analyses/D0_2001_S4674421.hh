// -*- C++ -*-
#ifndef RIVET_D0_2001_S4674421_HH
#define RIVET_D0_2001_S4674421_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {  


  /// @brief D0 Run I differential W/Z boson cross-section analysis
  /// @author Lars Sonnenschein
  class D0_2001_S4674421 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor.
    D0_2001_S4674421();    
    
    /// Factory method
    static Analysis* create() { 
      return new D0_2001_S4674421(); 
    }
    //@}
    
    
    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "4674421";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Run I differential W/Z boson cross-section analysis";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2001";
    }
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret += "arXiv:hep-ex/0107012";
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
    
    /// Analysis used ratio of mW/mZ 
    const double _mwmz;
    
    /// Ratio of \f$ BR(W->e,nu) \f$ used in the publication analysis
    const double _brwenu;
    
    /// Ratio of \f$ \text{BR}( Z \to e^+ e^-) \f$ used in the publication analysis
    const double _brzee;
    
    /// Invariant mass cuts for Z boson candidate (75 GeV < mZ < 105 GeV)
    const double _mZmin, _mZmax;


    // Event counters for cross section normalizations
    double _eventsFilledW;
    double _eventsFilledZ;
    
    //@{
    /// Histograms
    AIDA::IHistogram1D* _h_dsigdpt_w;
    AIDA::IHistogram1D* _h_dsigdpt_z;
    AIDA::IHistogram1D* _h_dsigdpt_scaled_z;
   //@}    

  };

}

#endif
