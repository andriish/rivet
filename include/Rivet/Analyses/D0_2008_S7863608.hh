// -*- C++ -*-
#ifndef RIVET_D0_2008_S7863608_HH
#define RIVET_D0_2008_S7863608_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/IsolationTools.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Measurement differential Z/gamma* + jet +X cross sections
  /// @author Gavin Hesketh
  class D0_2008_S7863608 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    D0_2008_S7863608();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7863608();
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "7863608";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Measurement of differential Z/gamma* + jet + X cross sections";
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
      ret += "Gavin Hesketh <gavin.hesketh@fnal.gov>";
      ret += "Frank Siegert <frank.siegert@durham.ac.uk>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Cross sections as a function of pT and rapidity of the boson "
         << "and pT and rapidity of the leading jet "
         << "in ppbar collisions at sqrt{s} = 1.96 TeV, based on "
         << "an integrated luminosity of 1.0 fb^-1.";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Tevatron Run 2 conditions:\n" << endl << endl
         << "* ppbar -> mu+ mu- + jets at 1960 GeV\n"
         << "* Needs mass cut on lepton pair to avoid photon singularity: min. range $65 < m_Z < 115$";
      return os.str();
    }
    string status() const {
      return "VALIDATED";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0808.1296 [hep-ex]");
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
    AIDA::IHistogram1D * _h_jet_pT_cross_section;
    AIDA::IHistogram1D * _h_jet_y_cross_section;
    AIDA::IHistogram1D * _h_Z_pT_cross_section;
    AIDA::IHistogram1D * _h_Z_y_cross_section;
    AIDA::IHistogram1D * _h_total_cross_section;
    //@}

  };


}

#endif
