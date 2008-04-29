// -*- C++ -*-
#ifndef RIVET_D0_2001_S4674421_HH
#define RIVET_D0_2001_S4674421_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/WZandh.hh"

namespace Rivet {  


  /// Implementation of D0 Run I, differential W/Z boson cross-section analysis.
  /// @author Lars Sonnenschein
  class D0_2001_S4674421 : public Analysis {

  public:

    /// Constructor.
    ///  - @c _mwmz = ratio of \f$ mW/mZ \f$ used in the publication analysis
    ///  - @c _brwenu = ratio of \f$ BR(W->e,nu) \f$ used in the publication analysis
    ///  - @c _brzee = ratio of \f$ BR(Z->ee) \f$ used in the publication analysis
    D0_2001_S4674421()
      : _mwmz(0.8820), _brwenu(0.1073), _brzee(0.033632)
    { 
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
      addProjection(*new WZandh(), "WZ");
    }    
    
    
    /// Factory method
    static Analysis* create() { 
      return new D0_2001_S4674421(); 
    }
    
    
    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "4674421";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Run I differential W/Z boson cross-section analysis";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2001";
    }
    /// Journal, and preprint references.
    vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("arXiv:hep-ex/0107012");
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
    
    /// analysis used ratio of mW/mZ 
    const double _mwmz;
    
    /// brwenu = ratio of \f$ BR(W->e,nu) \f$ used in the publication analysis
    const double _brwenu;
    
    /// brzee = ratio of \f$ BR(Z->ee) \f$ used in the publication analysis
    const double _brzee;
    
    // Event counters for cross section normalizations
    double _eventsFilledW;
    double _eventsFilledZ;
    
    //@{
    /// Histograms
    AIDA::IHistogram1D* _h_dsigdpt_w;
    AIDA::IHistogram1D* _h_dsigdpt_z;
    AIDA::IHistogram1D* _h_dsigdpt_wz_rat;
    //@}    

    /// Hide copy assignment operator.
    D0_2001_S4674421& operator=(const D0_2001_S4674421&);
  };

}

#endif
