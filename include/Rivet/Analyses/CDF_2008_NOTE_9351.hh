// -*- C++ -*-
#ifndef RIVET_CDF_2008_NOTE_9351_HH
#define RIVET_CDF_2008_NOTE_9351_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* @brief CDF Run II underlying event in Drell-Yan
   * @author Hendrik Hoeth
   * 
   * Measurement of the underlying event in Drell-Yan Z/gamma->e+e-
   * and Z/gamma->mu+mu- events. The reconstructed Z defines the
   * phi orientation. A Z mass window cut is applied.
   * 
   * 
   * @par Run conditions
   * 
   * @arg \f$ \sqrt{s} = \f$ 1960 GeV
   * @arg produce Drell-Yan events
   * @arg Set particles with c*tau > 10 mm stable
   * @arg Z decay mode: Z -> e+e- and Z -> mu+mu-
   * @arg gamma decay mode: gamma -> e+e- and gamma -> mu+mu-
   * @arg minimum invariant mass of the fermion pair coming from the Z/gamma: 70 GeV
   * 
   */ 
  class CDF_2008_NOTE_9351 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.5 \f$ GeV.
    CDF_2008_NOTE_9351();

    /// Factory method
    static Analysis* create() {
      return new CDF_2008_NOTE_9351();
    }
    //@}


  public:

    /// @name Publication metadata
    //@{

    /// Analysis name
    string name() const {
      return "CDF_2008_NOTE_9351";
    }
    /// SPIRES key (IRN)
    string spiresId() const {
      return "NONE";
    }
    /// A short description of the analysis.
    string summary() const {
      return "CDF Run 2 underlying event in Drell-Yan";
    }
    /// Full description of the analysis, for the manual
    string description() const {
      ostringstream os;
      os << "Deepak Kar's and Rick Field's measurement of the underlying event "
         << "in Drell-Yan events. Z -> ee and Z -> $\\mu\\mu$ events are selected "
         << "using a Z mass window cut between 70 and 110 GeV. ``Toward'', "
         << "``away'' and ``transverse'' regions are defined in the same way as "
         << "in the original (2001) CDF underlying event analysis. The reconstructed "
         << "Z defines the $\\phi$ direction of the toward region. The leptons are "
         << "ignored after the Z has been reconstructed. Thus the region most "
         << "sensitive to the underlying event is the toward region (the recoil jet "
         << "is boosted into the away region).";
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
     return "2008";
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
      os << "* Tevatron Run 2: ppbar collisions at 1960 GeV.\n"
         << "* Drell-Yan events with $Z/\\gamma* -> e e$ and $Z/\\gamma* -> \\mu\\mu$.\n"
         << "* A mass cut $m_{ll} > 70~\\text{GeV}$ can be applied on generator level.\n"
         << "* Particles with $c \\tau > 10~\\text{mm}$ should be set stable.";
      return os.str();
    }
    string status() const {
      return "VALIDATED";
    }
    /// Journal or preprint references
    vector<string> references() const {
      vector<string> ret;
      ret += "CDF public note 9351";
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

    AIDA::IProfile1D *_hist_tnchg;
    AIDA::IProfile1D *_hist_pnchg;
    AIDA::IProfile1D *_hist_pmaxnchg;
    AIDA::IProfile1D *_hist_pminnchg;
    AIDA::IProfile1D *_hist_pdifnchg;
    AIDA::IProfile1D *_hist_anchg;
    AIDA::IProfile1D *_hist_tcptsum;
    AIDA::IProfile1D *_hist_pcptsum;
    AIDA::IProfile1D *_hist_pmaxcptsum;
    AIDA::IProfile1D *_hist_pmincptsum;
    AIDA::IProfile1D *_hist_pdifcptsum;
    AIDA::IProfile1D *_hist_acptsum;
    AIDA::IProfile1D *_hist_tcptave;
    AIDA::IProfile1D *_hist_pcptave;
    AIDA::IProfile1D *_hist_acptave;
    AIDA::IProfile1D *_hist_tcptmax;
    AIDA::IProfile1D *_hist_pcptmax;
    AIDA::IProfile1D *_hist_acptmax;
    AIDA::IProfile1D *_hist_zptvsnchg;
    AIDA::IProfile1D *_hist_cptavevsnchg;
    AIDA::IProfile1D *_hist_cptavevsnchgsmallzpt;

  };


}

#endif
