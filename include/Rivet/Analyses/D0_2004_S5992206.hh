// -*- C++ -*-
#ifndef RIVET_D0_2004_S5992206_HH
#define RIVET_D0_2004_S5992206_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"

namespace Rivet {  


  /// @brief D0 Run II jet analysis
  /// @author Lars Sonnenschein
  class D0_2004_S5992206 : public Analysis {

  public:

    /// Constructor.
    D0_2004_S5992206() {
      setBeams(PROTON, ANTIPROTON);
      const FinalState fs(-3.0, 3.0);
      addProjection(fs, "FS");
      addProjection(D0ILConeJets(fs), "Jets");
      addProjection(TotalVisibleMomentum(fs), "CalMET");
      addProjection(PVertex(), "PV");

      // Veto neutrinos, and muons with pT above 1.0 GeV
      VetoedFinalState vfs(fs);
      vfs
        .addVetoPairId(NU_E)
        .addVetoPairId(NU_MU)
        .addVetoPairId(NU_TAU)
        .addVetoDetail(MUON, 1.0, MAXDOUBLE);
      addProjection(vfs, "VFS");
    }


    /// Factory method
    static Analysis* create() { 
      return new D0_2004_S5992206(); 
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "5992206";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Run II jet analysis";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2004";
    }
    /// Journal, and preprint references.
    vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("arXiv:hep-ex/0409040");
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

    /// Hide the assignment operator
    D0_2004_S5992206& operator=(const D0_2004_S5992206& x);

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
