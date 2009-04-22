// -*- C++ -*-
#ifndef RIVET_MC_TVT1960_ZJETS_HH
#define RIVET_MC_TVT1960_ZJETS_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Monte Carlo validation observables for Z[e+ e-] + jets production at Tevatron Run II
  /// @author Frank Siegert
  class MC_TVT1960_ZJETS : public Analysis {

  public:

    /// Default constructor.
    MC_TVT1960_ZJETS();


    /// Factory method 
    static Analysis* create() {
      return new MC_TVT1960_ZJETS();
    }


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "NONE";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Monte Carlo validation observables for Z[e+ e-] + jets production at Tevatron Run II";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "MC";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 2";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "NONE";
    }
    /// Overwrite virtual name method for special analysis scheme
    string name() const {
      return "MC_TVT1960_ZJETS";
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
      os << "Available observables\n" << endl << endl
         << "* Z mass\n" << endl
         << "* pT of jet 1-4\n" << endl
         << "* jet multiplicity\n" << endl
         << "* Delta eta (Z, jet1)\n" << endl
         << "* Delta R (jet2, jet3)\n" << endl
         << "* Differential jet rates 0->1, 1->2, 2->3, 3->4\n" << endl
         << "* Integrated 0-4 jet rates\n" << endl;
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Tevatron Run 2 conditions: " << endl << endl
         << "* ppbar -> e+ e- + jets at 1960 GeV. "
         << "* Needs mass cut on lepton pair to avoid photon singularity: "
         << "min. range $66 < m_ee < 116$ GeV" << endl;
      return os.str();
    }
    string status() const {
      return "NOT TO BE VALIDATED";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
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
    AIDA::IHistogram1D * _h_Z_mass;
    AIDA::IHistogram1D * _h_jet1_pT;
    AIDA::IHistogram1D * _h_jet2_pT;
    AIDA::IHistogram1D * _h_jet3_pT;
    AIDA::IHistogram1D * _h_jet4_pT;
    AIDA::IHistogram1D * _h_jet20_multi_exclusive;
    AIDA::IHistogram1D * _h_jet20_multi_inclusive;
    AIDA::IDataPointSet * _h_jet20_multi_ratio;
    AIDA::IHistogram1D * _h_jet10_multi_exclusive;
    AIDA::IHistogram1D * _h_jet10_multi_inclusive;
    AIDA::IDataPointSet * _h_jet10_multi_ratio;
    AIDA::IHistogram1D * _h_deta_Z_jet1;
    AIDA::IHistogram1D * _h_dR_jet2_jet3;
    AIDA::IHistogram1D * _h_log10_d[4];
    AIDA::IDataPointSet * _h_log10_R[5];
    //@}

  };

}

#endif
