// -*- C++ -*-
#ifndef RIVET_CDF_1994_S2952106_HH
#define RIVET_CDF_1994_S2952106_HH

#include "Rivet/Analysis.hh"

namespace Rivet {  

  /* @brief CDF Run I color coherence analysis
   * @author Lars Sonnenschein
   */
  class CDF_1994_S2952106 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_1994_S2952106();

    /// Factory method
    static Analysis* create() { 
      return new CDF_1994_S2952106(); 
    }
    //@}


    /// @name Publication metadata
    //@{

    /// SPIRES ID code.
    string spiresId() const {
      return "2952106";
    }

    /// A short description of the analysis.
    string summary() const {
      return "CDF Run I color coherence analysis.";
    }

    /// Full description of the analysis, to appear in the manual.
    string description() const {
      ostringstream os;
      os << "Analysis based on the CDF Run I color coherence analysis described in "
         << "PRD50,9,5562 (1994). Events with >= 3 jets are selected and Et "
         << "distributions of the three highest-pT jets are obtained. "
         << "$\\Delta{R}$ between 2nd and 3rd leading jets in pT and pseudorapidity "
         << "of the 3rd jet are plotted. $\\alpha = \\d{\\eta}/\\d{\\phi}$ is plotted, where "
         << "$\\d{\\eta}$ is the pseudorapidity difference between the 2nd and 3rd jets "
         << "and $\\d{\\phi}$ their azimuthal angle difference. Since the data has not been "
         << "detector-corrected, a bin by bin correction is applied, based on the "
         << "distributions with ideal and CDF simulation as given in the publication."
         << "NB. eta in [-4.2, 4.2] to match CDF CAL acceptance.";
      return os.str();
    }

    /// Event type required by this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Energy: sqrt(s) = 1800 GeV\n"
         << "* Event type: generic QCD events\n"
         << "* Cut on primary vertex $z$ position: z(PV) < 60 cm\n"
         << "* pTmin: leading jet = 100 GeV; and 3rd jet = 10 GeV\n"
         << "* Max. pseudorapidity range of 2nd and 3rd jets: $|eta| < 0.7$\n"
         << "* Azimuthal angle requirement: $\\Delta{\\phi} < \\pi/18$ (transverse back-to-backness)\n"
         << "* MET cut requirement: $mET/\\sqrt{\\text{Scalar }E_T} < 6.0 GeV$";
      return os.str();
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }

    /// Collider on which the experiment ran
    string collider() const {
      return "Tevatron Run 1";
    }

    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1994";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Lars Sonnenschein <Lars.Sonnenschein@cern.ch>";
      return rtn;
    }

    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Rev.D50,5562,1994";
      ret += "doi:10.1103/PhysRevD.50.5562";
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

    /// Counter for the number of events analysed
    double _eventsTried;

    /// Counter for the number of  3jet events passed
    double _events3jPassed;


    /// @name Analysis cuts
    //@{
    ///Cut on primary vertex z-position (z(PV) < 60 cm)
    const double _pvzmax;

    /// Min \f$ p_T \f$ of the leading and 3rd leading jets.
    //@{
    const double _leadJetPt;
    const double _3rdJetPt;
    //@}

    /// Max pseudorapidity range of 2nd and 3rd leading jets.
    const double _etamax;

    /// Delta phi (azimuthal angle) requirement (transverse back to back'ness).
    const double _phimin;

    /// MET over sqrt(scalar \f$ E_T \f$) cut requirement.
    const double _metsetmax;
    //@}


  private:

    /// @name Histogram collections
    //@{
    // AIDA::IHistogram2D* _histHvsDphi;
    // AIDA::IHistogram2D* _histRvsAlpha;
    AIDA::IHistogram1D* _histJet1Et;
    AIDA::IHistogram1D* _histJet2Et;
    AIDA::IHistogram1D* _histR23;
    AIDA::IHistogram1D* _histJet3eta;
    AIDA::IHistogram1D* _histAlpha;
    // AIDA::IHistogram1D* _histAlphaMCvsDat;
    AIDA::IHistogram1D* _histAlpaIdeal;
    AIDA::IHistogram1D* _histAlpaCDF;
    AIDA::IHistogram1D* _histR23Ideal;
    AIDA::IHistogram1D* _histR23CDF;
    AIDA::IHistogram1D* _histJet3etaIdeal;
    AIDA::IHistogram1D* _histJet3etaCDF;
    //@}

  };


}

#endif
