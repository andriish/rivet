// -*- C++ -*-
#ifndef RIVET_D0_2004_S5992206_HH
#define RIVET_D0_2004_S5992206_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/RivetAIDA.fhh"



namespace Rivet {  

  /// Analysis based on the D0 Run II jet analysis described in hep-ex/0409040.
  /// @author Lars Sonnenschein
  class D0_2004_S5992206 : public Analysis {

  public:

    /// Default constructor.
    D0_2004_S5992206()
      // NB. eta in [-3,3] cut specified via FinalState constructor
      : _fsproj(-3.0, 3.0), _vfsproj(_fsproj), 
	_conejetsproj(_fsproj), _calmetproj(_fsproj), _vertexproj()
    { 

      setBeams(PROTON, ANTIPROTON);


      // Add particle/antiparticle vetoing: 12=nu_e, 14=nu_mu, 16=nu_tau
      _vfsproj
        .addVetoPairId(12)
        .addVetoPairId(14)
        .addVetoPairId(16);
      
      // Veto muons (PDG code = 13) with pT above 1.0 GeV
      _vfsproj.addVetoDetail(13, 1.0, numeric_limits<double>::max());


      addProjection(_fsproj);
      addProjection(_vfsproj);
      addProjection(_conejetsproj);
      addProjection(_calmetproj);
      addProjection(_vertexproj);
      
    }

    /// Factory method
    static Analysis* create() { return new D0_2004_S5992206(); }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "5992206";
    }
    /// Get a description of the analysis.
    // string getDescription() const {
    //   return "";
    // }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2004";
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// The final state projector used by this analysis.
    FinalState _fsproj;

    ///The vetoed final state projector needed by the jet algorithm
    VetoedFinalState _vfsproj; 

    /// The D0ILConeJets projector used by this analysis.
    D0ILConeJets _conejetsproj;

    /// The Calorimeter Missing Et projector
    TotalVisibleMomentum _calmetproj;

    /// The Primary Vertex projector
    PVertex _vertexproj;

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
