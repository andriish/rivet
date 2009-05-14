// -*- C++ -*-
#ifndef RIVET_H1_2000_S4129130_HH
#define RIVET_H1_2000_S4129130_HH

#include "Rivet/Analysis.hh" 

namespace Rivet {


  /// @brief H1 energy flow and charged particle spectra
  /// @author Peter Richardson
  /// Based on the HZtool analysis
  class H1_2000_S4129130 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    H1_2000_S4129130();

    /// Factory method.
    static Analysis* create() { 
      return new H1_2000_S4129130(); 
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "4129130";
    }
    /// A short description of the analysis.
    string summary() const {
      return "H1 energy flow in DIS";
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Measurements of transverse energy flow for neutral current deep-inelastic "
         << "scattering events produced in positron-proton collisions at HERA. The kinematic "
         << "range covers squared momentum transfers $Q^2$ from 3.2 to 2200 GeV^2; the Bjorken "
         << "scaling variable $x$ from 8x10^-5 to 0.11 and the hadronic mass $W$ from "
         << "66 to 233 GeV. The transverse energy flow is measured in the hadronic centre of "
         << "mass frame and is studied as a function of $Q^2$, $x$, $W$ and pseudorapidity. "
         << "The behaviour of the mean transverse energy in the central pseudorapidity region "
         << "and an interval corresponding to the photon fragmentation region are analysed as "
         << "a function of Q^2 and W."
         << "\n\n"
         << "This analysis is useful for exploring the effect of photon PDFs and for tuning "
         << "models of parton evolution and treatment of fragmentation and the proton remnant " 
         << "in DIS.";
      return os.str();
    }
    /// Event type required by this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Event type: e+ p deep inelastic scattering\n"
         << "* Energy: p @ 820 GeV, e+ @ 27.5 GeV -> sqrt(s) = 300 GeV";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "H1";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "HERA";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2000";
    }
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Peter Richardson <peter.richardson@durham.ac.uk>";
      return rtn;
    }
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Eur.Phys.J.C12:595-607,2000";
      ret += "doi:10.1007/s100520000287";
      ret += "arXiv:hep-ex/9907027v1";
      return ret;
    }
    /// Validation status
    string status() const {
      return "VALIDATED";
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /**
     *  Polar angle with right direction of the beam
     */
    inline double beamAngle(const FourVector& v, const bool & order) {
      double thel = v.polarAngle()/degree;
      if(thel<0.) thel+=180.;
      if(!order) thel = 180.-thel;
      return thel;
    }

    /// @name Histograms
    //@{
    vector<AIDA::IHistogram1D *> _histETLowQa;
    vector<AIDA::IHistogram1D *> _histETHighQa;
    vector<AIDA::IHistogram1D *> _histETLowQb;
    vector<AIDA::IHistogram1D *> _histETHighQb;
    AIDA::IProfile1D * _histAverETCentral;
    AIDA::IProfile1D * _histAverETFrag;
    //@}

    /// @name storage of weights for normalisation
    //@{
    vector<double> _weightETLowQa;
    vector<double> _weightETHighQa;
    vector<double> _weightETLowQb;
    vector<double> _weightETHighQb;
    //@}
  };

}

#endif
