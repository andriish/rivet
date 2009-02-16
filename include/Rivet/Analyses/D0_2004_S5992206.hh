// -*- C++ -*-
#ifndef RIVET_D0_2004_S5992206_HH
#define RIVET_D0_2004_S5992206_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"

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
    /// A short description of the analysis.
    string spiresId() const {
      return "5992206";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Run II jet azimuthal decorrelation analysis";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2004";
    }
    /// Journal, and preprint references.
    vector<string> references() const {
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
