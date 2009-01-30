// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7782535_HH
#define RIVET_CDF_2008_S7782535_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// Implementation of CDF RunII b-jet shape paper
  /// @todo Test with Pythia
  class CDF_2008_S7782535 : public Analysis {

  public:

    /// Constructor.

    /// jet cuts: |eta| <= 0.7

    /// Don't attempt to model the following cuts :
    ///   Missing ET significance
    ///   veto on additional vertices
    ///   Zvtx < 50

    CDF_2008_S7782535()
      : _Rjet(0.7) , _NpTbins(4)
    { 
      setBeams(PROTON, ANTIPROTON);

      const FinalState fs(-3.6, 3.6);
      addProjection(fs, "FS");
      // Veto (anti)neutrinos, and muons with pT above 1.0 GeV
      VetoedFinalState vfs(fs);
      vfs
        .addVetoPairId(NU_E)
        .addVetoPairId(NU_MU)
        .addVetoPairId(NU_TAU)
        .addVetoDetail(MUON, 1.0*GeV, MAXDOUBLE);
      addProjection(vfs, "VFS");
      addProjection(FastJets(vfs, FastJets::CDFMIDPOINT, 0.7), "Jets");
      addProjection(JetShape(vfs, _jetaxes, 0.0, 0.7, 0.1, 0.3, ENERGY), "JetShape");
    }

  public:

    /// Factory method
    static Analysis* create() { 
      return new CDF_2008_S7782535(); 
    }


    /// @name Publication metadata
    //@{
    /// Get SPIRES ID code.
    string spiresId() const {
      return "7782535";
    }
    /// Get a description of the analysis.
    string description() const {
      return "CDF Run II b-jet shape paper";
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
    virtual vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0806.1699");
      return ret;
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}


  private:

    /// @name Analysis cuts

    /// @name Analysis data
    //@{
    vector<FourMomentum> _jetaxes;

    double _Rjet;
    vector<double> _pTbins;
    int _NpTbins;
    //@{
    /// Histograms
    AIDA::IProfile1D* _Psi_pT[4];
    AIDA::IDataPointSet* _OneMinusPsi_vs_pT;
    //@}

  };

}

#endif
