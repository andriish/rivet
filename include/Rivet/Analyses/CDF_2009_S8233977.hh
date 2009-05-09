// -*- C++ -*-
#ifndef RIVET_CDF_2009_S8233977_HH
#define RIVET_CDF_2009_S8233977_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* @brief CDF Run II min-bias cross-section
   * @author Hendrik Hoeth
   * 
   * Measurement of \f$ \langle p_T \rangle \f$ vs. \f$ n_\text{ch} \f$,
   * the track \f$ p_T \f$ distribution, and the \f$ \sum E_T \f$ distribution.
   * Particles are selected within |eta|<1 and with pT>0.4 GeV.
   * There is no pT cut for the \f$ \sum E_T \f$ measurement.
   * 
   * @par Run conditions
   * 
   * @arg \f$ \sqrt{s} = \f$ 1960 GeV
   * @arg Run with generic QCD events.
   * @arg Set particles with c*tau > 10 mm stable
   * 
   */ 
  class CDF_2009_S8233977 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.4 \f$ GeV.
    CDF_2009_S8233977();

    /// Factory method
    static Analysis* create() {
      return new CDF_2009_S8233977();
    }
    //@}


  public:

    /// @name Publication metadata
    //@{

    /// SPIRES key (IRN)
    string spiresId() const {
      return "8233977";
    }
    /// A short description of the analysis.
    string summary() const {
      return "CDF Run 2 min bias cross-section analysis";
    }
    /// Full description of the analysis, for the manual
    string description() const {
      ostringstream os;
      os << "Niccolo' Moggi's minbias analysis. Minimum bias events are used to "
         << "measure the average track $p_T$ vs charged multiplicity, a track $p_T$ "
         << "distribution and an inclusive $\\sum E_T$ distribution. "
         << "\n\n "
         << "WARNING: Only average track $p_T$ vs charged multiplicity is validated!";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
     return "CDF";
    }
    /// Collider on which the experiment was based
    string collider() const {
     return "Tevatron Run 2";
    }
    /// When published according to SPIRES
    string year() const {
     return "2009";
    }
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Hendrik Hoeth <hendrik.hoeth@cern.ch>";
      return ret;
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Tevatron Run 2: ppbar QCD interactions at 1960 GeV.\n"
         << "* Particles with $c \\tau > {}$10 mm should be set stable.";
      return os.str();
    }

    string status() const {
      return "PARTIALLY VALIDATED";
    }
    /// No journal or preprint references: this is a demo.
    vector<string> references() const {
      vector<string> ret;
      ret += "CDF public note 9337";
      ret += "hep-ex/0904.1098";
      return ret;
    }
    //@}


  public:

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    AIDA::IProfile1D *_hist_pt_vs_multiplicity;
    AIDA::IHistogram1D *_hist_pt;
    AIDA::IHistogram1D *_hist_sumEt;

  };


}

#endif
