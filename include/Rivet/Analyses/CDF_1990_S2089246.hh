// -*- C++ -*-
#ifndef RIVET_CDF_1990_S2089246_HH
#define RIVET_CDF_1990_S2089246_HH

#include "Rivet/Analysis.hh"

namespace Rivet {  

  /* @brief CDF pseudorapidity analysis
   * @author Andy Buckley
   */
  class CDF_1990_S2089246 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_1990_S2089246();

    /// Factory method
    static Analysis* create() { 
      return new CDF_1990_S2089246(); 
    }
    //@}


    /// @name Publication metadata
    //@{

    /// SPIRES ID code.
    string spiresId() const {
      return "2089246";
    }

    /// A short description of the analysis.
    string summary() const {
      return "CDF pseudorapidity distributions at 630 and 1800 GeV";
    }

    /// Full description of the analysis, to appear in the manual.
    string description() const {
      ostringstream os;
      os << "Pseudorapidity distributions based on the CDF 630 "
         << "and 1800 GeV runs from 1987. All data is detector "
         << "corrected. The data confirms the UA5 measurement of "
         << "a $\\d{N}/\\d{\\eta}$ rise with energy faster than "
         << "$\\ln{\\sqrt{s}}$, and as such this analysis is important "
         << "for constraining the energy evolution of minimum bias "
         << "and underlying event characteristics in MC simulations.";
      return os.str();
    }

    /// Event type required by this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Energy: sqrt(s) = 630 and 1800 GeV\n"
         << "* Event type: generic QCD events\n"
         << "* $|\\eta| < 3.5$";
      return os.str();
    }

    /// Validation status
    string status() const {
      return "UNVALIDATED";
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }

    /// Collider on which the experiment ran
    string collider() const {
      return "Tevatron Run 0";
    }

    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1990";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Andy Buckley <andy.buckley@cern.ch>";
      return rtn;
    }

    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Rev.D41:2330,1990";
      ret += "doi:10.1103/PhysRevD.41.2330";
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

    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _hist_eta630;
    AIDA::IHistogram1D* _hist_eta1800;
    //@}

  };


}

#endif
