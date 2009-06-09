// -*- C++ -*-
#ifndef RIVET_CDF_2008_LEADINGJETS_HH
#define RIVET_CDF_2008_LEADINGJETS_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* CDF Run II underlying event in leading jet events
   * @author Hendrik Hoeth
   * 
   * Rick Field's measurement of the underlying event in "leading jet" events.
   * The leading jet (CDF midpoint R=0.7) must be within |eta|<2 and defines
   * the "toward" phi direction. Particles are selected in |eta|<1. For the pT
   * related observables there is a pT>0.5 GeV cut. For sum(ET) there is no pT cut.
   * 
   * 
   * @par Run conditions
   * 
   * @arg \f$ \sqrt{s} = \f$ 1960 GeV
   * @arg Run with generic QCD events.
   * @arg Set particles with c*tau > 10 mm stable
   * @arg Several \f$ p_\perp^\text{min} \f$ cutoffs are probably required to fill the profile histograms:
   *   @arg \f$ p_\perp^\text{min} = \f$ 0 (min bias), 10, 20, 50, 100, 150 GeV
   *   @arg The corresponding merging points are at \f$ p_T = \f$ 0, 30, 50, 80, 130, 180 GeV
   * 
   */ 
  class CDF_2008_LEADINGJETS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    CDF_2008_LEADINGJETS();

    /// Factory method
    static Analysis* create() {
      return new CDF_2008_LEADINGJETS();
    }
    //@}


  public:

    /// @name Publication metadata
    //@{

    /// Analysis name
    string name() const {
      return "CDF_2008_LEADINGJETS";
    }
    /// SPIRES key (IRN)
    string spiresId() const {
      return "NONE";
    }
    /// A short description of the analysis.
    string summary() const {
      return "CDF Run 2 underlying event in leading jet events";
    }
    /// Full description of the analysis, for the manual
    string description() const {
      ostringstream os;
      os << "Rick Field's measurement of the underlying event in leading jet "
         << "events. If the leading jet of the event is within $|\\eta| < 2$, "
	 << "the event is accepted and ``toward'', ``away'' and ``transverse'' "
	 << "regions are defined in the same way as in the original (2001) CDF "
	 << "underlying event analysis. The leading jet defines the $\\phi$ "
	 << "direction of the toward region. The transverse regions are most "
	 << "sensitive to the underlying event.";
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
      os << "Tevatron Run 2: ppbar QCD interactions at 1960 GeV. "
         << "Particles with $c \\tau > {}$10 mm should be set stable. "
	 << "Several $p_\\perp^\\text{min}$ cutoffs are probably required to "
	 << "fill the profile histograms: "
	 << " * $p_\\perp^\\text{min} = {}$ 0 (min bias), 10, 20, 50, 100, 150 GeV "
         << " * The corresponding merging points are at $p_T = $ 0, 30, 50, 80, 130, 180 GeV";
      return os.str();
    }
    /// Validation status
    string status() const {
      return "VALIDATED";
    }
    /// No journal or preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret += "";
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

    AIDA::IProfile1D *_hist_pnchg;
    AIDA::IProfile1D *_hist_pmaxnchg;
    AIDA::IProfile1D *_hist_pminnchg;
    AIDA::IProfile1D *_hist_pdifnchg;
    AIDA::IProfile1D *_hist_pcptsum;
    AIDA::IProfile1D *_hist_pmaxcptsum;
    AIDA::IProfile1D *_hist_pmincptsum;
    AIDA::IProfile1D *_hist_pdifcptsum;
    AIDA::IProfile1D *_hist_pcptave;

  };


}

#endif
