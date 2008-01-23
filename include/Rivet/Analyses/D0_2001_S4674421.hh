// -*- C++ -*-
#ifndef RIVET_D0_2001_S4674421_HH
#define RIVET_D0_2001_S4674421_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/WZandh.hh"
//#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {  

  /// Implementation of D0 Run I, differential W/Z boson cross 
  /// section analyisis, publication hep-ex/0107012 
  /// _mwmz = ratio of \f$ mW/mZ \f$ used in the publication analysis
  /// _brwenu = ratio of \f$ BR(W->e,nu) \f$ used in the publication analysis
  /// _brzee = ratio of \f$ BR(Z->ee) \f$ used in the publication analysis
  class D0_2001_S4674421 : public Analysis {

  public:

    /// Constructor.
    D0_2001_S4674421()
      : _WZproj(), 
        _mwmz(0.8820), _brwenu(0.1073), _brzee(0.033632)
    { 
      
      setBeams(PROTON, ANTIPROTON);
      addProjection(_WZproj);
      setNeedsCrossSection(true);
    }    
    
    
    /// Factory method
    static Analysis* create() { return new D0_2001_S4674421(); }
    
    
    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "4674421";
    }
    /// Get a description of the analysis.
    string getDescription() const {
       return "";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2001";
    }
    //@}
    
    
    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}
    
  private:
    
    /// The final state projector used by this analysis.
    WZandh _WZproj;
    
    /// analysis used ratio of mW/mZ 
    const double _mwmz;
    
    /// brwenu = ratio of \f$ BR(W->e,nu) \f$ used in the publication analysis
    const double _brwenu;
    
    /// brzee = ratio of \f$ BR(Z->ee) \f$ used in the publication analysis
    const double _brzee;
    
    /// Hide copy assignment operator.
    D0_2001_S4674421& operator=(const D0_2001_S4674421&);
    
    // Event counters for cross section normalizations
    double _eventsFilledW;
    double _eventsFilledZ;
    
    //@{
    /// Histograms
    AIDA::IHistogram1D* _h_dsigdpt_w;
    AIDA::IHistogram1D* _h_dsigdpt_z;
    AIDA::IHistogram1D* _h_dsigdpt_wz_rat;
    //@}    
    
  };

}

#endif
