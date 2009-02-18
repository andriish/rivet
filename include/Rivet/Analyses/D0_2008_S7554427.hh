// -*- C++ -*-
#ifndef RIVET_D0_2008_S7554427_HH
#define RIVET_D0_2008_S7554427_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"

namespace Rivet {


  /// @brief Measurement of D0 Run II Z pT differential cross-section shape
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7554427 : public Analysis {

  public:

    /// Default constructor.
    D0_2008_S7554427();


    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7554427();
    }


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "7554427";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Z/gamma* + X cross-section shape, differential in pT(Z)";
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
      ret += "Andy Buckley <andy.buckley@durham.ac.uk>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Boson pT in Z/gamma* [e+ e-] + X  events " << endl
         << "======================================== " << endl
         << "Cross sections as a function of pT of the vector boson "
         << "inclusive and in forward region (|y|>2, pT<30GeV)"
         << "in ppbar collisions at sqrt{s} = 1.96 TeV, based on "
         << "an integrated luminosity of 0.98 fb^-1." << endl;
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Tevatron Run 2 conditions: " << endl << endl
         << "* ppbar -> e+ e- + jets at 1960 GeV. "
         << "* Needs mass cut on lepton pair to avoid photon singularity: min. range 40 < m_ee < 200" << endl;
      return os.str();
    }
    string status() const {
      return "VALIDATED";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("hep-ex/0712.0803");
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
    AIDA::IHistogram1D * _h_ZpT;
    AIDA::IHistogram1D * _h_forward_ZpT;
    //@}

  };


}

#endif
