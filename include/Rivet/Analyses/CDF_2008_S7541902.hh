// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7541902_HH
#define RIVET_CDF_2008_S7541902_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief CDF jet pT and multiplicity distributions in W + jets events
  ///
  /// This CDF analysis provides jet pT distributions for 4 jet multiplicity bins
  /// as well as the jet multiplicity distribution in W + jets events.
  /// e-Print: arXiv:0711.4044 [hep-ex]
  class CDF_2008_S7541902 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2008_S7541902();

    /// Factory method.
    static Analysis* create() { 
      return new CDF_2008_S7541902(); 
    }

    //@}


    /// @name Publication metadata
    //@{
    /// Get the SPIRES ID
    string spiresId() const {
      return "7541902";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Jet pT distributions for 4 jet multiplicity bins as well as the jet multiplicity distribution in W + jets events.";
    }
    
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os <<  "Measurement of the cross section for W boson production in association with jets in $p\\bar{p}$ " 
	 <<  "collisions at $\\sqrt{s}=1.96$ TeV. The analysis uses 320 pb$^{-1}$ of data collected with the "
	 <<  "CDF II detector. "
	 <<  "W bosons are identified in their $e\\nu$ decay channel and "
	 <<  "jets are reconstructed using an $R < 0.4$ cone algorithm. "
	 <<  "For each $W + \\geq$ n-jet sample (where n = 1--4) a measurement of "
	 <<  "d$\\sigma({p}\\bar{p} \\rightarrow W + \\geq$ n jet)/d$E_T(n^{th}$-jet) $\\times$ BR($W \\rightarrow{e}\\nu$) "
	 <<  "is made, where d$E_T(n^{th}$-jet) is the Et of the n$^{th}$-highest Et jet above 20 GeV. "
	 <<  "A measurement of the total cross section, " 
	 <<  "$\\sigma(p\\bar{p} \\rightarrow W + \\geq$ n-jet) $\\times$ BR($W \\rightarrow{e}\\nu)$ with "
	 <<  "$E_T(n^{th}-jet) > 25$ GeV is also made. " 
	 <<  "Both measurements are made for jets with $|\\eta| < 2$ and for a limited reigon of "
	 <<  "the $W \\rightarrow{e}\\nu$ decay phase space: $|\\eta^{e}| < 1.1$, $p_{T}^{e} > 20$ GeV, "
	 <<  "$p_{T}^{\\nu} > 30$ GeV and $M_{T} > 20$ GeV. "
	 <<  "The cross sections are corrected for all detector effects and can be directly compared "
	 <<  "to particle level $W$ + jet(s) predictions. "
	 <<  "These measurements can be used to test and tune QCD predictions for the number of jets "
	 <<   "in and kinematics of $W$ + jets events.";
      return os.str();
    }

    /// Characteristics of events to be processed by this analysis
    string runInfo() const {
      ostringstream os;
      os << "Requires the process $p\\bar{p} \\rightarrow {W} "
         <<  "\\rightarrow{e}\\nu$, "
         << "additional hard jets will also have to be included to get a good description. " 
         << "The LO process in Herwig is set with IPROC=1451.";
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
    /// Collider on which the experiment ran.
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
      rtn += "Ben Cooper <b.d.cooper@qmul.ac.uk>";
      rtn += "Emily Nurse <nurse@hep.ucl.ac.uk>";
      return rtn;
    }
    
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0711.4044 [hep-ex]");
      ret.push_back("Phys.Rev.D77:011108,2008");
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

    /// @name Cuts 
    //@{
    /// Cut on the electron ET:
    double _electronETCut;
    /// Cut on the electron ETA:
    double _electronETACut;   
    /// Cut on the missing ET
    double _eTmissCut;
    /// Cut on the transverse mass squared
    double _mTCut;
    /// Cut on the jet ET for differential cross sections
    double _jetEtCutA;
    /// Cut on the jet ET for jet multiplicity
    double _jetEtCutB;
    /// Cut on the jet ETA
    double _jetETA;
    //@}    

    double _xpoint;

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histJetEt[4];
    AIDA::IHistogram1D* _histJetMultNorm;
    AIDA::IDataPointSet* _histJetMultRatio[4];
    AIDA::IHistogram1D* _histJetMult[4];
    //@}




  };


}

#endif
