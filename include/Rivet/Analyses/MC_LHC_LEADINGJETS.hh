// -*- C++ -*-
#ifndef RIVET_MC_LHC_LEADINGJETS_HH
#define RIVET_MC_LHC_LEADINGJETS_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* Underlying event in leading jet, extended to the LHC
   * @author Andy Buckley
   */ 
  class MC_LHC_LEADINGJETS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_LHC_LEADINGJETS();

    /// Factory method
    static Analysis* create() {
      return new MC_LHC_LEADINGJETS();
    }
    //@}


  public:

    /// @name Publication metadata
    //@{

    /// Analysis name
    string name() const {
      return "MC_LHC_LEADINGJETS";
    }
    /// SPIRES key (IRN)
    string spiresId() const {
      return "NONE";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Underlying event in leading jet events, extended to LHC";
    }
    /// Full description of the analysis, for the manual
    string description() const {
      ostringstream os;
      os << "Rick Field's measurement of the underlying event in leading jet "
         << "events, extended to the LHC. As usual, the leading jet of the "
         << "defines an azimuthal toward/transverse/away decomposition, "
         << "in this case the event is accepted within $|\\eta| < 2$, as "
         << "in the CDF 2008 version of the analysis. Since this isn't the "
         << "Tevatron, I've chosen to use $k_\\perp$ rather than midpoint "
         << "jets.";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
     return "NONE";
    }
    /// Collider on which the experiment was based
    string collider() const {
     return "LHC";
    }
    /// When published according to SPIRES
    string year() const {
     return "NONE";
    }
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Andy Buckley <andy.buckley@cern.ch>";
      return ret;
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "LHC: pp QCD interactions at 0.9, 10 or 14 TeV."
         << "Particles with $c \\tau > {}$10 mm should be set stable."
         << "Several $p_\\perp^\\text{min}$ cutoffs are probably required to "
         << "fill the profile histograms. ";
      //<< " * $p_\\perp^\\text{min} = {}$ 0 (min bias), 10, 20, 50, 100, 150 GeV "
      //<< " * The corresponding merging points are at $p_T = $ 0, 30, 50, 80, 130, 180 GeV";
      return os.str();
    }

    string status() const {
      return "NOT TO BE VALIDATED";
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
    AIDA::IProfile1D *_hist_pcptsum;
    AIDA::IProfile1D *_hist_pmaxcptsum;
    AIDA::IProfile1D *_hist_pmincptsum;
    AIDA::IProfile1D *_hist_pcptave;

  };


}

#endif
