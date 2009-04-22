// -*- C++ -*-
#ifndef RIVET_D0_2004_S5992206_HH
#define RIVET_D0_2004_S5992206_HH

#include "Rivet/Analysis.hh"

namespace Rivet {  


  /* @brief D0 Run II jet analysis
   * @author Lars Sonnenschein
   * 
   * Measurement of angular correlations in di-jet events.
   * 
   * 
   * @par Run conditions
   * 
   * @arg \f$ \sqrt{s} = \f$ 1960 GeV
   * @arg Run with generic QCD events.
   * @arg Several \f$ p_\perp^\text{min} \f$ cutoffs are probably required to fill the histograms:
   *   @arg \f$ p_\perp^\text{min} = \f$ 50, 75, 100, 150 GeV for the four pT ranges respecively
   * 
   */ 
  class D0_2004_S5992206 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor.
    D0_2004_S5992206();

    /// Factory method
    static Analysis* create() { 
      return new D0_2004_S5992206(); 
    }

    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "5992206";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Run II jet azimuthal decorrelation analysis";
    }
    /// Full description of the analysis, for the manual
    string description() const {
      ostringstream os;
      os << "Correlations in the azimuthal angle between the two largest $p_T$ jets "
         << "have been measured using the D0 detector in ppbar collisions at 1960 GeV. "
	 << "The analysis is based on an inclusive dijet event sample in the central "
	 << "rapidity region. The correlations are determined for four different $p_T$ "
	 << "intervals.";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "D0";
    }
    /// Collider on which the experiment was based
    string collider() const {
      return "Tevatron Run 2";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2004";
    }
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Lars Sonnenschein <lars.sonnenschein@cern.ch>";
      return ret;
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Tevatron Run 2: ppbar QCD interactions at 1960 GeV.";
      return os.str();
    }
    string status() const {
      return "VALIDATED";
    }
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys. Rev. Lett., 94, 221801 (2005)";
      ret += "arXiv:hep-ex/0409040";
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
    AIDA::IHistogram1D* _histJetAzimuth_pTmax75_100;
    AIDA::IHistogram1D* _histJetAzimuth_pTmax100_130;
    AIDA::IHistogram1D* _histJetAzimuth_pTmax130_180;
    AIDA::IHistogram1D* _histJetAzimuth_pTmax180_;
    //@}

  };

}

#endif
