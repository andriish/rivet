// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7541902_HH
#define RIVET_CDF_2008_S7541902_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/SVertex.hh"

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
    CDF_2008_S7541902()
      : _electronETCut(20.0 *GeV), _electronETACut(1.1),
        _eTmissCut(30.0 *GeV), _mTCut(20.0 *GeV),
        _jetEtCutA(20.0 *GeV),  _jetEtCutB(25.0 *GeV), _jetETA(2.0),
        _xpoint(1960.)
    {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);

      // Basic FS
      FinalState fs(-3.6, 3.6);
      addProjection(fs, "FS");

      // Create a final state with any e-nu pair with invariant mass 65 -> 95 GeV and ET > 20 (W decay products)
      std::vector<std::pair<long,long> > vids;
      vids.push_back(make_pair(ELECTRON, NU_EBAR));
      vids.push_back(make_pair(POSITRON, NU_E));
      FinalState fs2(-3.6, 3.6, 20*GeV);
      InvMassFinalState invfs(fs2, vids, 65*GeV, 95*GeV);
      addProjection(invfs, "INVFS");

      // Make a final state without the W decay products for jet clustering
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(invfs);
      addProjection(vfs, "VFS");
      addProjection(FastJets(vfs, FastJets::CDFJETCLU, 0.4), "Jets");
    }


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
      os <<  "Measurement of the cross section for W boson production in association with jets in $p\bar{p}$" 
	 <<  "collisions at $\sqrt{s}=1.96$ TeV. The analysis uses 320 pb$^{-1}$ of data collected with the"
	 <<  "CDF II detector."
	 <<  "W bosons are identified in their $e\nu$ decay channel and"
	 <<  "jets are reconstructed using an $R < 0.4$ cone algorithm."
	 <<  "For each $W + \geq$ n-jet sample (where n = 1--4) a measurement of"
	 <<  "d$\sigma({p}\bar{p} \rightarrow W + \geq$ n jet)/d$E_T(n^{th}$-jet) $\times$ BR($W \rightarrow{e}\nu$)"
	 <<  "is made, where d$E_T(n^{th}$-jet) is the Et of the n$^{th}$-highest Et jet above 20 GeV."
	 <<  "A measurement of the total cross section," 
	 <<  "$\sigma(p\bar{p} \rightarrow W + \geq$ n-jet) $\times$ BR($W \rightarrow{e}\nu)$ with"
	 <<  "$E_T(n^{th}-jet) > 25$ GeV is also made." 
	 <<  "Both measurements are made for jets with $|\eta| < 2$ and for a limited reigon of"
	 <<  "the $W \rightarrow{e}\nu$ decay phase space: $|\eta^{e}| < 1.1$, $p_{T}^{e} > 20$ GeV,"
	 <<  "$p_{T}^{\nu} > 30$ GeV and $M_{T} > 20$ GeV."
	 <<  "The cross sections are corrected for all detector effects and can be directly compared"
	 <<  "to particle level $W$ + jet(s) predictions."
	 <<  "These measurements can be used to test and tune QCD predictions for the number of jets"                             <<   "in and kinematics of $W$ + jets events.";
      return os.str();
    }
  /// Characteristics of events to be processed by this analysis
    string runInfo() const {
      ostringstream os;
      os << "Requires the process ppbar->W->enu,"
         << "additional hard jets will also have to be included to get a good description." 
	 << "The LO process in Herwig is set with IPROC=1451";
      return os.str();
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
