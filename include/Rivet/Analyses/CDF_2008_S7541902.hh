// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7541902_HH
#define RIVET_CDF_2008_S7541902_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/SVertex.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  
  /// This CDF analysis provides jet pT distributions for 4 jet multplicity bins
  /// as well as the jet multiplicity distribution in W + jets events.
  /// e-Print: arXiv:0711.4044 [hep-ex]
  class CDF_2008_S7541902 : public Analysis {

  public:

    /// Default constructor
    /// _fsproj     : final state particles
    /// _vfsproj    : _fsproj without neutrinos, muons with \f$ p_T \f$ above 1.0 GeV and W decay products for jet clustering
    /// _conejetsproj   : jetClu clusters with R = 0.4, running on _vfsproj
    CDF_2008_S7541902()
      : _fsproj(), _vfsproj(_fsproj), _conejetsproj(_vfsproj, FastJets::CDFJETCLU, 0.4)
    {

      setBeams(PROTON, ANTIPROTON);

      // Add particle/antiparticle vetoing for jet clustering: 12=nu_e, 14=nu_mu, 16=nu_tau
      _vfsproj
        .addVetoPairId(12)
	.addVetoPairId(14)
	.addVetoPairId(16);
      // Veto muons (PDG code = 13) with \f$ p_T \f$ above 1.0 GeV
      _vfsproj.addVetoDetail(13, 1.0, numeric_limits<double>::max());
      // Veto the W decay products 
      _vfsproj.addDecayProductsVeto(24);
      _vfsproj.addDecayProductsVeto(-24);
      addProjection(_fsproj);
      addProjection(_vfsproj);
      addProjection(_conejetsproj);
      setNeedsCrossSection(true);
    }

    /// Factory method.
    static Analysis* create() { 
      return new CDF_2008_S7541902(); 
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "7541902";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2008";
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// The final state projector.
    FinalState _fsproj;
    /// The visible final state projector.
    VetoedFinalState _vfsproj;
    /// The FastJets projector used by this analysis.
    FastJets _conejetsproj;

  private:
    
    /// Hide the assignment operator
    CDF_2008_S7541902& operator=(const CDF_2008_S7541902&);

    //@{
    /// Histograms
    AIDA::IHistogram1D* _histJetMult;
    AIDA::IHistogram1D* _histJetEt[4];
    //@}
  };

}

#endif
