// -*- C++ -*-
#ifndef RIVET_D0_2009_S8320160_HH
#define RIVET_D0_2009_S8320160_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class D0_2009_S8320160 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    D0_2009_S8320160();

    /// Factory method 
    static Analysis* create() {
      return new D0_2009_S8320160();
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis. 
    string spiresId() const {
      return "8320160";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Dijet angular distributions";
    }    
    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "D0";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 2";
    }
    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2009";
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
      os << "Dijet angular distributions in different bins of dijet mass "
         << "from 0.25 TeV to above 1.1 TeV "
         << "in ppbar collisions at $\\sqrt{s} = 1.96$ TeV, based on "
         << "an integrated luminosity of $0.7 \\text{fb}^{-1}$.";
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
      ret.push_back("arXiv:0906.4819 [hep-ex]");
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
    BinnedHistogram<double> _h_chi_dijet;
    //@}
    
  };


}

#endif
