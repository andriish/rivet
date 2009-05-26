// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7540469_HH
#define RIVET_CDF_2008_S7540469_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement differential Z/gamma* + jet + X cross sections
  /// @author Frank Siegert
  class CDF_2008_S7540469 : public Analysis {

  public:

    /// Default constructor.
    CDF_2008_S7540469();


    /// Factory method 
    static Analysis* create() {
      return new CDF_2008_S7540469();
    }


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "7540469";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Measurement of differential Z/gamma* + jet + X cross sections";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "CDF";
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
      ret += "Frank Siegert <frank.siegert@durham.ac.uk>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Cross sections as a function of jet transverse momentum in 1 and "
         << "2 jet events, and jet multiplicity in ppbar collisions at sqrt(s) "
         << "= 1.96 TeV, based on an integrated luminosity of 1.7 fb^-1. The "
         << "measurements cover the rapidity region $|y_\\text{jet}| < 2.1$ "
         << "and the transverse momentum range $pT_\\text{jet} > 30~\\text{GeV}/c$." << endl;
      return os.str();
    }

    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Tevatron Run 2 conditions: ppbar -> e+ e- + jets at 1960 GeV.\n"
         << "* Needs mass cut on lepton pair to avoid photon singularity: min. range $66 < m_ee < 116$" << endl;
      return os.str();
    }

    string status() const {
      return "VALIDATED";
    }

    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Rev.Lett.100:102001,2008";
      //ret += "doi:TODO";
      ret += "arXiv:0711.3717 [hep-ex]";
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
    AIDA::IHistogram1D * _h_jet_multiplicity;
    AIDA::IHistogram1D * _h_jet_pT_cross_section_incl_1jet;
    AIDA::IHistogram1D * _h_jet_pT_cross_section_incl_2jet;
    //@}

  };

}

#endif
