// -*- C++ -*-
#ifndef RIVET_ALEPH_1991_S2435284_HH
#define RIVET_ALEPH_1991_S2435284_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement of ALEPH LEP1 charged multiplicity
  /// @author Andy Buckley
  class ALEPH_1991_S2435284 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor.
    ALEPH_1991_S2435284();

    /// Factory method.
    static Analysis* create() { 
      return new ALEPH_1991_S2435284(); 
    }
    //@}


  public:

    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "2435284";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Hadronic Z decay charged multiplicity measurement";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "ALEPH";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "LEP 1";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1991";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Andy Buckley <andy.buckley@durham.ac.uk>";
      return ret;
    }

    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "The charged particle multiplicity distribution of hadronic Z decays, as measured "
         << "on the peak of the Z resonance using the ALEPH detector at LEP. The unfolding procedure "
         << "was model independent, and the distribution was found to have a mean of 20.85+-0.24, "
         << "Comparison with lower energy data supports the KNO scaling hypothesis. " 
         << "The shape of the multiplicity distribution is well described by a log-normal "
         << "distribution, as predicted from a cascading model for multi-particle production. ";
      return os.str();
    }

    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Hadronic Z decay events generated on the Z pole (sqrt(s) = 91.2 GeV)";
      return os.str();
    }

    string status() const {
      return "VALIDATED";
    }

    /// Journal, and preprint references
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys. Lett. B, 273, 181 (1991)";
      return ret;
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event & event);
    virtual void finalize();
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histChTot;
    //@}

  };

}


#endif
