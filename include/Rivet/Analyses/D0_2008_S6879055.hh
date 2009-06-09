// -*- C++ -*-
#ifndef RIVET_D0_2008_S6879055_HH
#define RIVET_D0_2008_S6879055_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement of the ratio sigma(Z/gamma* + n jets)/sigma(Z/gamma*)
  class D0_2008_S6879055 : public Analysis {

  public:

    /// Default constructor.
     D0_2008_S6879055();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S6879055();
    }

    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "6879055";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Measurement of the ratio sigma(Z/gamma* + n jets)/sigma(Z/gamma*)";
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
      return "2008";
    }
    /// Names & emails of analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Giulio Lenzi";
      ret += "Frank Siegert <frank.siegert@durham.ac.uk>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Cross sections as a function of pT of the three leading jets "
         << "and n-jet cross section ratios in ppbar collisions at sqrt{s} "
         << "= 1.96 TeV, based on an integrated luminosity of 0.4 fb^-1.";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Tevatron Run 2 conditions:\n"
         << "* ppbar -> e+ e- + jets at 1960 GeV\n"
         << "* Needs mass cut on lepton pair to avoid photon singularity: min. range $75 < m_ee < 105$";
      return os.str();
    }
    /// Validation status
    string status() const {
      return "VALIDATED";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("hep-ex/0608052");
      return ret;
    }

    //@}
    /// @name Analysis methods
    //@{ 
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _crossSectionRatio;
    AIDA::IHistogram1D * _pTjet1;
    AIDA::IHistogram1D * _pTjet2;
    AIDA::IHistogram1D * _pTjet3;
    //@}

  };

}

#endif
