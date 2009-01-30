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
        _eTmissCut(30.0 *GeV), _mT2Cut(200.0 * GeV * GeV),
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
    /// Get a description of the analysis.
    string description() const {
      return "Jet pT distributions for 4 jet multiplicity bins as well as the jet multiplicity distribution in W + jets events.";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2008";
    }
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0711.4044 [hep-ex]");
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
    double _mT2Cut;
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
