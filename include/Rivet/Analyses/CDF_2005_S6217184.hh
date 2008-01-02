// -*- C++ -*-
#ifndef RIVET_CDF_2005_S6217184_HH
#define RIVET_CDF_2005_S6217184_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/JetShape.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// Implementation of CDF RunII jet shape paper hep-ex/0505013 
  class CDF_2005_S6217184 : public Analysis {

  public:

    /// Constructor.
    /// \f$ \eta \in [-2,2] \f$ cut used on final state.
    /// Jet shape \f$ r_\text{min} = 0.0 \f$, \f$ r_\text{max} = 0.7 \f$, 
    /// interval = 0.1, r1minPsi = 0.3.
    CDF_2005_S6217184()
      : _fsproj(-2.0, 2.0), _vfsproj(_fsproj), _jetsproj(_fsproj), 
        _calmetproj(_fsproj), _vertexproj(),
        _jetshapeproj(_vfsproj, _jetaxes, 0.0, 0.7, 0.1, 0.3, ENERGY)
    { 
      setBeams(PROTON, ANTIPROTON);

      // Add particle/antiparticle vetoing: 12=nu_e, 14=nu_mu, 16=nu_tau
      /// @todo Use ParticleName enum for clarity.
      _vfsproj
        .addVetoPairId(12)
        .addVetoPairId(14)
        .addVetoPairId(16);
      
      // Veto muons (PDG code = 13) with pT above 1.0 GeV
      /// @todo Use ParticleName enum for clarity.
      _vfsproj.addVetoDetail(13, 1.0, numeric_limits<double>::max());

      addProjection(_fsproj);
      addProjection(_vfsproj);
      addProjection(_jetsproj);
      addProjection(_calmetproj);
      addProjection(_vertexproj);
      addProjection(_jetshapeproj);

      _Rjet = 0.7;

      /// @todo Ahem - replace this.
      _pTbins.resize(19);
      _pTbins[0]  =  37.0;
      _pTbins[1]  =  45.0;
      _pTbins[2]  =  55.0;
      _pTbins[3]  =  63.0;
      _pTbins[4]  =  73.0;
      _pTbins[5]  =  84.0;
      _pTbins[6]  =  97.0;
      _pTbins[7]  = 112.0;
      _pTbins[8]  = 128.0;
      _pTbins[9]  = 148.0;
      _pTbins[10] = 166.0;
      _pTbins[11] = 186.0;
      _pTbins[12] = 208.0;
      _pTbins[13] = 229.0;
      _pTbins[14] = 250.0;
      _pTbins[15] = 277.0;
      _pTbins[16] = 304.0;
      _pTbins[17] = 340.0;
      _pTbins[18] = 380.0;
      for (size_t i = 0; i < 18; ++i) _ShapeWeights[i] = 0.0;
    }
      

  public:

    /// Factory method
    static Analysis* create() { return new CDF_2005_S6217184(); }


    /// @name Publication metadata
    //@{
    /// Get SPIRES ID code.
    string getSpiresId() const {
      return "6217184";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "CDF Run II jet shape paper";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2005";
    }
    /// Journal, and preprint references.
    virtual vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("hep-ex/0505013");
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

    /// The final state projector used by this analysis.
    FinalState _fsproj;
    
    ///The vetoed final state projector needed by the jet algorithm.
    VetoedFinalState _vfsproj; 
    
    /// The FastJets projector is used by this analysis.
    FastJets _jetsproj;
    
    /// The calorimeter missing \f$ E_T \f$ projector.
    TotalVisibleMomentum _calmetproj;

    /// The primary vertex projector.
    PVertex _vertexproj;

    /// The jet shape projector.
    JetShape _jetshapeproj;
    

  private:

    vector<FourMomentum> _jetaxes;

    double _Rjet;

    double _ShapeWeights[18];

    /// pT bins to be distiguished during analysis
    vector<double> _pTbins;

    /// Hide the assignment operator
    CDF_2005_S6217184& operator=(const CDF_2005_S6217184&);

    //@{
    /// Histograms
    AIDA::IProfile1D* _profhistRho_pT[18];
    AIDA::IProfile1D* _profhistPsi_pT[18];
    AIDA::IProfile1D* _profhistPsi;
    //@}

  };

}

#endif
