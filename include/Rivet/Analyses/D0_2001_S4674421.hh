// -*- C++ -*-
#ifndef RIVET_D0_2001_S4674421_HH
#define RIVET_D0_2001_S4674421_HH

#include "Rivet/Analysis.hh"

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
      return "Tevatron Run I differential W/Z boson cross-section analysis";
    }

    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Measurement of differential W/Z boson cross section and ratio in p pbar "
         << "collisions at center-of-mass energy sqrt(s) = 1.8 TeV. The data cover "
         << "electrons and neutrinos in a pseudo-rapidity range of -2.5 to 2.5. "
        //   << "\n\n"
         << "";
      return os.str();
    }

    /// Event type required by this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Energy: sqrt(s) = 1800 GeV\n"
         << "* Event type: W/Z events with decays to first generation leptons"
         << "";
      return os.str();
    }

    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "D0";
    }

    /// Collider on which the experiment ran
    string collider() const {
      return "Tevatron Run 1";
    }

    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2001";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Lars Sonnenschein <Lars.Sonnenschein@cern.ch>";
      return rtn;
    }

    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Lett.B517:299-308,2001";
      ret += "doi:10.1016/S0370-2693(01)01020-6";
      ret += "arXiv:hep-ex/0107012v2";
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
