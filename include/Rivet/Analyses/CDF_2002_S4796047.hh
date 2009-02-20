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


    /// SPIRES ID code.
    string spiresId() const {
      return "4796047";
    }

    /// A short description of the analysis.
    string summary() const {
      return "CDF Run 1 charged multiplicity measurement";
    }

    /// Full description of the analysis, to appear in the manual.
    string description() const {
      ostringstream os;
      os << "A study of pp̅ collisions at sqrt(s) = 1800 and 630 GeV collected using a minimum bias "
         << "trigger in which the data set is divided into two classes corresponding to `soft' and " 
         << "`hard' interactions. For each subsample, the analysis includes measurements of the "
         << "multiplicity, transverse momentum (pT) spectra, and the average pT and event-by-event "
         << "pT dispersion as a function of multiplicity. A comparison of results shows distinct "
         << "differences in the behavior of the two samples as a function of the center of mass "
         << "energy. The properties of the soft sample are invariant as a function of c.m. energy.";
      return os.str();
    }

    /// Event type required by this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Energy: sqrt(s) = 630 and 1800 GeV\n"
         << "* Event type: generic QCD events\n"
         << "* TODO: MORE?";
      return os.str();
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }

    /// Collider on which the experiment ran
    string collider() const {
      return "Tevatron Run 1";
    }

    /// When published (according to SPIRES). 
    string year() const {
      return "2002";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Hendrik Hoeth <hendrik.hoeth@cern.ch>";
      return rtn;
    }

    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Rev.D65:072005,2002";
      ret += "doi:10.1103/PhysRevD.65.072005 ";
      return ret;
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

  };


}

#endif
