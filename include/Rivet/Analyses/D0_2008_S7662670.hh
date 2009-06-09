// -*- C++ -*-
#ifndef RIVET_D0_2008_S7662670_HH
#define RIVET_D0_2008_S7662670_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement of D0 differential jet cross sections
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7662670 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_2008_S7662670();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7662670();
    }
    //@}


    /// @name Publication metadata
    //@{

    /// A short description of the analysis. 
    string spiresId() const {
      return "7662670";
    }

    /// A short description of the analysis.
    string summary() const {
      return "Measurement of D0 Run II differential jet cross sections";
    }

    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Measurement of the inclusive jet cross section in p pbar collisions at "
         << "center-of-mass energy sqrt(s) = 1.96 TeV. The data cover jet transverse "
         << "momenta from 50--600 GeV and jet rapidities in the range -2.4 to 2.4. "
        //   << "\n\n"
         << "";
      return os.str();
    }

    /// Event type required by this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Energy: sqrt(s) = 1960 GeV\n"
         << "* Event type: QCD events\n"
         << "* pTmin cut may be necessary: lowest jet pT bin is at 50 GeV";
      return os.str();
    }

    /// Validation status
    string status() const {
      return "VALIDATED";
    }

    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "D0";
    }

    /// Collider on which the experiment ran
    string collider() const {
      return "Tevatron Run 2";
    }

    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2008";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Andy Buckley <andy.buckley@durham.ac.uk>";
      rtn += "Gavin Hesketh <gavin.hesketh@cern.ch>";
      return rtn;
    }

    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Rev.Lett.101:062001,2008";
      ret += "doi:10.1103/PhysRevLett.101.062001";
      ret += "arXiv:0802.2400v3 [hep-ex]";
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
    AIDA::IHistogram1D* _h_dsigdptdy_y00_04;
    AIDA::IHistogram1D* _h_dsigdptdy_y04_08;
    AIDA::IHistogram1D* _h_dsigdptdy_y08_12;
    AIDA::IHistogram1D* _h_dsigdptdy_y12_16;
    AIDA::IHistogram1D* _h_dsigdptdy_y16_20;
    AIDA::IHistogram1D* _h_dsigdptdy_y20_24;
    //@}

  };

}

#endif
