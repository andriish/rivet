// -*- C++ -*-
#ifndef RIVET_CDF_2005_S6217184_HH
#define RIVET_CDF_2005_S6217184_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /* @brief CDF Run II jet shape analysis
   * @author Lars Sonnenschein
   * @author Andy Buckley
   */	
  class CDF_2005_S6217184 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    // Constructor
    /// Cuts on final state are \f$ \eta \in [-2,2] \f$.
    /// Jet shape \f$ r_\text{min} = 0.0 \f$, \f$ r_\text{max} = 0.7 \f$, 
    /// interval = 0.1, r1minPsi = 0.3.
    CDF_2005_S6217184();

    /// Factory method
    static Analysis* create() { 
      return new CDF_2005_S6217184(); 
    }

    //@}


    /// @name Publication metadata
    //@{

    /// SPIRES ID code.
    string spiresId() const {
      return "6217184";
    }

    /// A short description of the analysis.
    string summary() const {
      return "CDF Run II jet shape analysis";
    }

    /// Full description of the analysis, to appear in the manual.
    string description() const {
      ostringstream os;
      os << "Measurement of jet shapes in inclusive jet production in p pbar collisions at"
         << "center-of-mass energy sqrt(s) = 1.96 TeV. The data cover jet transverse "
         << "momenta from 37--380 GeV and absolute jet rapidities in the range 0.1 to 0.7. ";
      return os.str();
    }

    /// Event type required by this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Energy: sqrt(s) = 1960 GeV\n"
         << "* Event type: generic QCD events.\n"
         << "* $\\eta \\in [-2,2]$ cut used on final state.\n"
         << "* Jet axes must have $|y| \\in [0.1, 0.7]$.\n"
         << "* Jet shape $r \\in [0.0, 0.7]$\n"
         << "* Jet pTmin in plots is 37 GeV/c: choose generator min pT somewhere well below this.";
      return os.str();
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }

    /// Collider on which the experiment ran
    string collider() const {
      return "Tevatron Run 2";
    }

    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2005";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Lars Sonnenschein <Lars.Sonnenschein@cern.ch>";
      rtn += "Andy Buckley <andy.buckley@cern.ch>";
      return rtn;
    }

    /// Journal, and preprint references.
    virtual vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Rev.D71:112002,2005";
      ret += "doi:10.1103/PhysRevD.71.112002";
      ret += "arXiv:hep-ex/0505013";
      return ret;
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}


  private:

    /// @name Analysis data
    //@{

    /// Vector of jet axes
    vector<FourMomentum> _jetaxes;

    /// \f$p_\perp\f$ bins to be distinguished during analysis
    vector<double> _pTbins;
    //@}


    /// @name Histograms
    //@{
    AIDA::IProfile1D* _profhistRho_pT[18];
    AIDA::IProfile1D* _profhistPsi_pT[18];
    AIDA::IProfile1D* _profhistPsi;
    //@}

  };


}

#endif
