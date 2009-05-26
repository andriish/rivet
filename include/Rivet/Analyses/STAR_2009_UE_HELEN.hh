// -*- C++ -*-
#ifndef RIVET_STAR_2009_UE_HELEN_HH
#define RIVET_STAR_2009_UE_HELEN_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* STAR underlying event
   * @author Hendrik Hoeth
   */ 
  class STAR_2009_UE_HELEN : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    STAR_2009_UE_HELEN();

    /// Factory method
    static Analysis* create() {
      return new STAR_2009_UE_HELEN();
    }
    //@}


  public:

    /// @name Publication metadata
    //@{

    /// Analysis name
    string name() const {
      return "STAR_2009_UE_HELEN";
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
      os << "";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
     return "STAR";
    }
    /// Collider on which the experiment was based
    string collider() const {
     return "(RHIC pp 200 GeV)";
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
      os << "* pp interactions at 200 GeV";
      return os.str();
    }

    string status() const {
      return "UNVALIDATED";
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
