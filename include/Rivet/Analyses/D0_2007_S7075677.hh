// -*- C++ -*-
#ifndef RIVET_D0_2007_S7075677_HH
#define RIVET_D0_2007_S7075677_HH

#include "Rivet/Analysis.hh"

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
    /// A short description of the analysis. 
    string spiresId() const {
      return "7075677";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Z/gamma* + X cross-section shape, differential in y(Z)";
    }
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "D0";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 2";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2007";
    }
    /// Names & emails of analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Andy Buckley <andy.buckley@durham.ac.uk>";
      ret += "Gavin Hesketh <ghesketh@fnal.gov>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Cross sections as a function of boson rapidity "
         << "ppbar collisions at sqrt{s} = 1.96 TeV, based on "
         << "an integrated luminosity of 0.4 fb^-1.";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Tevatron Run 2 conditions:\n"
         << "* ppbar -> e+ e- + jets at 1960 GeV.\n"
         << "* Needs mass cut on lepton pair to avoid photon singularity: min. range $71 < m_ee < 111$";
      return os.str();
    }
    string status() const {
      return "UNCLEAR: Photons in Z reconstruction?";
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
