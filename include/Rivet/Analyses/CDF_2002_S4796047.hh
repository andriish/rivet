// -*- C++ -*-
#ifndef RIVET_CDF_2002_S4796047_HH
#define RIVET_CDF_2002_S4796047_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /*
   * @brief CDF Run I charged multiplicity measurement
   * @author Hendrik Hoeth
   * 
   * This analysis measures the charged multiplicity distribution
   * in minimum bias events at two different center-of-mass energies:
   * \f$ \sqrt{s} = \f$ 630 and 1800 GeV.
   * 
   * Particles with c*tau > 10 mm are considered stable, i.e. they
   * are reconstructed and their decay products removed. Selection
   * cuts are |eta|<1 and pT>0.4 GeV.
   * 
   * 
   * @par Run conditions
   * 
   * @arg Two different beam energies: \f$ \sqrt{s} = \$f 630 & 1800 GeV
   * @arg Run with generic QCD events.
   * @arg Set particles with c*tau > 10 mm stable
   * 
   */
  class CDF_2002_S4796047 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.4 \f$ GeV.
    CDF_2002_S4796047()
    { 
      setBeams(PROTON, ANTIPROTON);
      addProjection(Beam(), "Beams");
      const ChargedFinalState cfs(-1.0, 1.0, 0.4*GeV);
      addProjection(cfs, "FS");
    }


    /// Factory method
    static Analysis* create() {
      return new CDF_2002_S4796047();
    }
    //@}


  public:

    string getSpiresId() const {
      return "4796047";
    }
    string getDescription() const {
      return "CDF Run 1 charged multiplicity measurement";
    }
    string getExpt() const {
      return "CDF";
    }
    string getYear() const {
      return "2002";
    }

  public:

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    AIDA::IHistogram1D *_hist_multiplicity_630;
    AIDA::IHistogram1D *_hist_multiplicity_1800;

  private:

    /// Hide the assignment operator.
    CDF_2002_S4796047& operator=(const CDF_2002_S4796047&);

  };

}

#endif
