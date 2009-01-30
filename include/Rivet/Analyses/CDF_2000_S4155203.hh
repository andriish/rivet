// -*- C++ -*-
#ifndef RIVET_CDF_2000_S4155203_HH
#define RIVET_CDF_2000_S4155203_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {


  /*
   * @brief CDF Run I Z pT in Drell-Yan events
   * @author Hendrik Hoeth
   * 
   * Measurement of the Z pT distribution in Z -> e+e- events
   * at a center-of-mass energy of \f$ \sqrt{s} = \f$ 1800 GeV.
   * A Z mass window cut is applied.
   * 
   * 
   * @par Run conditions
   * 
   * @arg \f$ \sqrt{s} = \f$ 1800 GeV
   * @arg produce Drell-Yan events
   * @arg Z decay mode: Z -> e+e-
   * @arg gamma decay mode: gamma -> e+e-
   * @arg minimum invariant mass of the fermion pair coming from the Z/gamma: 66 GeV
   * 
   */ 
  class CDF_2000_S4155203 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.5 \f$ GeV.
    CDF_2000_S4155203()
    { 
      setBeams(PROTON, ANTIPROTON);
      const ChargedFinalState clfs(-4.2, 4.2, 15*GeV);
      addProjection(ChargedLeptons(clfs), "CL");
    }


    /// Factory method
    static Analysis* create() {
      return new CDF_2000_S4155203();
    }
    //@}


  public:

    string spiresId() const {
      return "4155203";
    }
    string description() const {
      return "CDF Run 1 Z pT measurement in Z->e+e- events";
    }
    string experiment() const {
      return "CDF";
    }
    string year() const {
      return "2000";
    }


  public:

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    AIDA::IHistogram1D *_hist_zpt;

  };


}

#endif
