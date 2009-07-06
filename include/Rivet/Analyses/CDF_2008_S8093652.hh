// -*- C++ -*-
#ifndef RIVET_CDF_2008_S8093652_HH
#define RIVET_CDF_2008_S8093652_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CDF_2008_S8093652 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    CDF_2008_S8093652();

    /// Factory method 
    static Analysis* create() {
      return new CDF_2008_S8093652();
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "8093652";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Dijet mass spectrum";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "CDF";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 2";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2008";
    }
    /// Names & emails of analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Frank Siegert <frank.siegert@durham.ac.uk>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Dijet mass spectrum  "
         << "from 0.2 TeV to 1.4 TeV "
         << "in ppbar collisions at $\\sqrt{s} = 1.96$ TeV, based on "
         << "an integrated luminosity of $1.13 \\text{fb}^{-1}$.";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Tevatron Run 2 conditions:\n" << endl << endl
         << "* ppbar -> jets at 1960 GeV";
      return os.str();
    }
    string status() const {
      return "UNVALIDATED";
    }
    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0812.4036 [hep-ex]");
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
    AIDA::IHistogram1D* _h_m_dijet;
    //@}
    
  };


}

#endif
