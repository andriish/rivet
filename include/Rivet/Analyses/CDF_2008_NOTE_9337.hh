// -*- C++ -*-
#ifndef RIVET_CDF_2008_NOTE_9337_HH
#define RIVET_CDF_2008_NOTE_9337_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

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
  class CDF_2008_NOTE_9337 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.4 \f$ GeV.
    CDF_2008_NOTE_9337()
    { 
      setBeams(PROTON, ANTIPROTON);
      const FinalState fs(-1.0, 1.0, 0.0*GeV);
      const ChargedFinalState cfs(-1.0, 1.0, 0.4*GeV);
      addProjection(fs, "FS");
      addProjection(cfs, "CFS");
      setNeedsCrossSection(true);
    }


    /// Factory method
    static Analysis* create() {
      return new CDF_2008_NOTE_9337();
    }
    //@}


  public:

    string name() const {
      return "CDF_2008_NOTE_9337";
    }
    string description() const {
      return "CDF Run 2 min bias cross-section analysis";
    }

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

  private:

    /// Hide the assignment operator.
    CDF_2008_NOTE_9337& operator=(const CDF_2008_NOTE_9337&);

  };

}

#endif
