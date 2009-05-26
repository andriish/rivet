// -*- C++ -*-
#ifndef RIVET_D0_2008_S7719523_HH
#define RIVET_D0_2008_S7719523_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement of isolated gamma + jet + X differential cross-sections
  /// Inclusive isolated gamma + jet cross-sections, differential in pT(gamma), for 
  /// various photon and jet rapidity bins.
  ///
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7719523 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
     D0_2008_S7719523();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7719523();
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "7719523";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Isolated gamma + jet cross-sections, differential in pT(gamma) for various y-bins";
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "The process p pbar -> photon + jet + X as studied by the D0 detector at "
         << "the Fermilab Tevatron collider at center-of-mass energy sqrt(s) = 1.96 TeV. "
         << "Photons are reconstructed in the central rapidity region $|y_\\gamma| < 1.0$ "
         << "with transverse momenta in the range 30--400 GeV, while jets are reconstructed "
         << "in either the central $|y_jet| < 0.8$ or forward $1.5 < |y_\\text{jet}| < 2.5$ "
         << "rapidity intervals with $pT_\\text{jet} > 15~\\text{GeV}$. The differential cross section "
         << "$\\mathrm{d}^3 \\sigma / \\mathrm{d}{pT_\\gamma} \\mathrm{d}{y_\\gamma} "
         << "\\mathrm{d}{y_jet}$ is measured as a function of $pT_\\gamma$ in four regions, "
         << "differing by the relative orientations of the photon and the jet."
        //<< "Ratios between the differential cross sections in each region are also presented."
         << "\n\n"
         << "MC predictions have trouble with simultaneously describing the measured "
         << "normalization and $pT_\\gamma$ dependence of the cross section in any of the "
         << "four measured regions.";
      return os.str();
    }
    /// Details of the type of events expected by this analysis
    string runInfo() const {
      ostringstream os;
      os << "* Produce only gamma + jet (q,qbar,g) hard processes\n\n"
         << "  * for Pythia 6: MSEL=10, and MSUB indices 14, 29 & 115 enabled\n\n"
         << "* Lowest bin edge at 30 GeV: kinematic pTmin cut may be required for "
         << "  good statistics.";
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
      ret += "Phys.Lett.B666:435-445,2008";
      ret += "doi:10.1016/j.physletb.2008.06.076";
      ret += "arXiv:0804.1107v2 [hep-ex]";
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
    AIDA::IHistogram1D* _h_central_same_cross_section;
    AIDA::IHistogram1D* _h_central_opp_cross_section;
    AIDA::IHistogram1D* _h_forward_same_cross_section;
    AIDA::IHistogram1D* _h_forward_opp_cross_section;
    //@}

  };

}

#endif
