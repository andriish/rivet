// -*- C++ -*-
#ifndef RIVET_STAR_2008_S7993412_HH
#define RIVET_STAR_2008_S7993412_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief di-hadron correlations in d-Au at 200 GeV
  class STAR_2008_S7993412 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    STAR_2008_S7993412();

    /// Factory method 
    static Analysis* create() {
      return new STAR_2008_S7993412();
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "7993412";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Di-hadron correlations in d-Au at 200 GeV";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "STAR";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "RHIC d-Au 200 GeV";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2008";
    }
    /// Names & emails of analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Christine Nattrass <christine.nattrass@yale.edu>";
      ret += "Hendrik Hoeth <hendrik.hoeth@cern.ch>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Correlation in $\\eta$ and $\\phi$ between the "
         << "charged hadron with the highest $p_T$ (\"trigger "
         << "particle\") and the other charged hadrons in the "
         << "event (\"associated particles\"). The data was "
         << "collected in d-Au collisions at 200 GeV. "
         << "Nevertheless, it is very proton-proton like and "
         << "can therefore be compared to pp Monte Carlo "
         << "(not for tuning, but for qualitative studies).";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "RHIC d-Au 200 GeV run conditions:\n" << endl << endl
         << "* d-Au at 200 GeV (use pp Monte Carlo! See description.)";
      return os.str();
    }
    string status() const {
      return "UNVALIDATED";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0809.5261 [nucl-ex]");
      return ret;
    }
    //@}


    /// @name Analysis methods
    //@{ 
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IProfile1D * _h_Y_jet_trigger;
    AIDA::IProfile1D * _h_Y_jet_associated;
    //@}

  };


}

#endif
