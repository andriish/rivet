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

    string summary() const {
      return "Z pT measurement in CDF Z -> e+ e- events";
    }

    string description() const {
      ostringstream os;
      os << "Measurement of transverse momentum and total cross section of e^+e^- pairs "
         << "in the Z-boson region of 66 < M_{ee} < 116 GeV/c^2 from pbar-p collisions "
         << "at sqrt(s) = 1.8 TeV, with the Tevatron CDF detector."
         << "\n\n"
         << "The Z pT, in a fully-factorised picture, is generated by the momentum balance "
         << "against initial state radiation (ISR) and the primordial/intrinsic pT of the "
         << "Z's parent partons in the incoming hadrons. The Z pT is important in generator "
         << "tuning to fix the interplay of ISR and multi-parton interactions (MPI) in"
         << "generating 'underlying event' activity.";
      return os.str();
    }

    string runInfo() const {
      ostringstream os;
      os << "* Event type: Z Drell-Yan with e+ e- decay mode only.";
      return os.str();
    }

    string experiment() const {
      return "CDF";
    }

    string collider() const {
      return "Tevatron Run 1";
    }

    string year() const {
      return "2000";
    }

    vector<string> references() const {
      vector<string> rtn;
      rtn += "Phys.Rev.Lett.84:845-850,2000";
      rtn += "arXiv:hep-ex/0001021";
      rtn += "doi:10.1103/PhysRevLett.84.845";
      return rtn;
    }

    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Hendrik Hoeth <hendrik.hoeth@cern.ch>";
      return rtn;
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
