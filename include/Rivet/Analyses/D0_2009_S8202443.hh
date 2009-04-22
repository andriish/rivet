// -*- C++ -*-
#ifndef RIVET_D0_2009_S8202443_HH
#define RIVET_D0_2009_S8202443_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  class D0_2009_S8202443 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    D0_2009_S8202443();

    /// Factory method 
    static Analysis* create() {
      return new D0_2009_S8202443();
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "8202443";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Z/gamma* + jet + X cross sections differential in pT(jet 1,2,3)";
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
      return "2009";
    }
    /// Names & emails of analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Frank Siegert <frank.siegert@durham.ac.uk>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Cross sections as a function of pT of the three leading jets "
         << "in $Z/\\gamma^{*} (\\to e^{+} e^{-})$ + jet + X production "
         << "in ppbar collisions at $\\sqrt{s} = 1.96$ TeV, based on "
         << "an integrated luminosity of $1.0 fb^{-1}$.";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Tevatron Run 2 conditions:\n" << endl << endl
         << "* ppbar -> e+ e- + jets at 1960 GeV\n"
         << "* Needs mass cut on lepton pair to avoid photon singularity: min. range $65 < m_Z < 115$";
      return os.str();
    }
    string status() const {
      return "NOT VALIDATED";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0903.1748 [hep-ex]");
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
    AIDA::IHistogram1D * _h_jet1_pT;
    AIDA::IHistogram1D * _h_jet2_pT;
    AIDA::IHistogram1D * _h_jet3_pT;
    //@}

  };


}

#endif
