// -*- C++ -*-
#ifndef RIVET_CDF_2008_NOTE_9351_HH
#define RIVET_CDF_2008_NOTE_9351_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

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
    CDF_2008_NOTE_9351()
    { 
      setBeams(PROTON, ANTIPROTON);
      const ChargedFinalState cfs(-1.0, 1.0, 0.5*GeV);
      const ChargedFinalState clfs(-1.0, 1.0, 20*GeV);
      addProjection(cfs, "FS");
      addProjection(ChargedLeptons(clfs), "CL");
    }


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
      os << "CDF Run 2 underlying event in Drell-Yan. "
         << "TODO: MORE!";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
     return "CDF Run 2";
    }
    /// Collider on which the experiment was based
    string collider() const {
     return "Tevatron";
    }
    /// When published according to SPIRES
    string year() const {
     return "2008";
    }
    /// No journal or preprint references: this is a demo.
    vector<string> references() const {
      vector<string> ret;
      ret += "CDF/PUB/CDF/PUBLIC/9351";
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
