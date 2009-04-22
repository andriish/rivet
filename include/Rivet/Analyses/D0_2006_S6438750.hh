// -*- C++ -*-
#ifndef RIVET_D0_2006_S6438750_HH
#define RIVET_D0_2006_S6438750_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Inclusive isolated photon cross-section, differential in \f$ p_\perp(gamma) \f$.
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2006_S6438750 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Default constructor.
     D0_2006_S6438750();

    /// Factory method 
    static Analysis* create() {
      return new D0_2006_S6438750();
    }

    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "6438750";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Inclusive isolated photon cross-section, differential in pT(gamma)";
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Measurement of differential cross section for inclusive production of isolated photons "
         << "in p pbar collisions at sqrt(s) = 1.96 TeV with the D0 detector at the Fermilab Tevatron "
         << "collider. The photons span transverse momenta 23--300 GeV and have pseudorapidity "
         << "$|eta| < 0.9$."
         << "\n\n"
         << "Isolated direct photons are probes of pQCD via the annihilation (q qbar -> gamma g) "
         << "and quark-gluon Compton scattering (q g -> gamma q) processes, the latter of which is "
         << "also sensitive to the gluon PDF. The initial state radiation / resummation formalisms "
         << "are sensitive to the resulting photon pT spectrum";
      return os.str();
    }
    /// Characteristics of events to be processed by this analysis
    string runInfo() const {
      ostringstream os;
      os << "* Requires gamma + jet (q,qbar,g) hard processes\n\n"
         << "  * for Pythia 6, MSEL=10 for with MSUB indices 14, 18, 29, 114, 115 enabled\n\n"
         << "* Lowest pT bin is at 23 GeV: a pT_min cut at ~10--15 GeV may be "
         << "  required to get good statistics.";
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
      return "2006";
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
      ret += "Phys.Lett.B639:151-158,2006, Erratum-ibid.B658:285-289,2008";
      ret += "doi:10.1016/j.physletb.2006.04.048";
      ret += "arXiv:hep-ex/0511054 (plus erratum)";
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
    AIDA::IHistogram1D* _h_pTgamma;
    //@}

  };

}

#endif
