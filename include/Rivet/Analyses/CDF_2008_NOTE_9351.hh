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

    string name() const {
      return "CDF_2008_NOTE_9351";
    }
    string description() const {
      return "CDF Run 2 underlying event in Drell-Yan";
    }

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

  private:

    /// Hide the assignment operator.
    CDF_2008_NOTE_9351& operator=(const CDF_2008_NOTE_9351&);

  };

}

#endif
