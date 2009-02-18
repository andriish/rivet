// -*- C++ -*-
#ifndef RIVET_D0_2008_S7837160_HH
#define RIVET_D0_2008_S7837160_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"

namespace Rivet {


  /// @brief Measurement of W charge asymmetry from D0 Run II
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7837160 : public Analysis {

  public:

    /// Default constructor.
    D0_2008_S7837160();


    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7837160();
    }


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "7837160";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Measurement of W charge asymmetry from D0 Run II";
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Measurement of the electron charge asymmetry in p pbar -> W + X -> e nu_e  + X "
         << "events at a center of mass energy of 1.96 TeV. The asymmetry is measured as " 
         << "a function of the electron transverse momentum and pseudorapidity in the "
         << "interval (-3.2, 3.2)."
         << "\n\n"
         << "This data is sensitive to proton parton distribution functions due to the "
         << "valence asymmetry in the incoming quarks which produce the W. Initial state "
         << "radiation should also affect the pT distribution.";
      return os.str();
    }
    /// Type of events required by this analysis
    string runInfo() const {
      ostringstream os;
      os << "* Event type: W production with decay to e nu_e only\n\n"
         << "  * for Pythia 6: MSEL = 12, MDME(206,1) = 1\n\n"
         << "* Energy: 1.96 TeV";
      return os.str();
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
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Andy Buckley <andy.buckley@durham.ac.uk>";
      rtn += "Gavin Hesketh <gavin.hesketh@cern.ch>";
      return rtn;
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Rev.Lett.101:211801,2008";
      ret += "doi:10.1103/PhysRevLett.101.211801";
      ret += "arXiv:0807.3367v1 [hep-ex]";
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
    /// @todo Move into histo manager
    AIDA::IHistogram1D *_h_dsigplus_deta_25_35, *_h_dsigminus_deta_25_35;
    AIDA::IHistogram1D *_h_dsigplus_deta_35, *_h_dsigminus_deta_35;
    AIDA::IHistogram1D *_h_dsigplus_deta_25, *_h_dsigminus_deta_25;
    //@}

  };


}

#endif
