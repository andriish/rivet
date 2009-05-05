// -*- C++ -*-
#ifndef RIVET_STAR_2006_S6870392_HH
#define RIVET_STAR_2006_S6870392_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief inclusive jet cross-section in pp at 200 GeV
  class STAR_2006_S6870392 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    STAR_2006_S6870392();

    /// Factory method 
    static Analysis* create() {
      return new STAR_2006_S6870392();
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "6870392";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Inclusive jet cross-section in pp at 200 GeV";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "STAR";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "RHIC pp 200 GeV";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2006";
    }
    /// Names & emails of analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Hendrik Hoeth <hendrik.hoeth@cern.ch>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Inclusive jet cross section as a function of $p_T$ "
         << "in $pp$ collisions at $\\sqrt{s} = 200$ GeV, measured "
         << "by the STAR experiment at RHIC.";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "RHIC pp 200 GeV run conditions:\n" << endl << endl
         << "* pp at 200 GeV";
      return os.str();
    }
    string status() const {
      return "VALIDATED";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("Phys. Rev. Lett. 97, 252001");
      ret.push_back("hep-ex/0608030");
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
    AIDA::IHistogram1D * _h_jet_pT_MB;
    AIDA::IHistogram1D * _h_jet_pT_HT;
    //@}

  };


}

#endif
