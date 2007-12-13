// -*- C++ -*-
#ifndef RIVET_CDF_2006_S6653332_HH
#define RIVET_CDF_2006_S6653332_HH

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

  
  /// This CDF analysis provides \f$ p_T \f$ and \f$ \eta \f$ distributions 
  /// of jets in Z + (b) jet production, before and after tagging.
  class CDF_2006_S6653332 : public Analysis {

  public:

    /// Default constructor
    /// Max. distance (eta,phi) between vertex vis. momentum and 
    /// jet to be probed = 0.7
    /// Tracker geometrical acceptance (eta < 2.0)
    /// Impact Parameter resolution = 34e-3mm, including beam spot
    /// cut on Decay Length Significance = 7.5
    /// Decay Length Significance resolution (assumed to be 34e-3mm)
    CDF_2006_S6653332()
      : _fsproj(-3.6, 3.6), _vfsproj(_fsproj), _jetsproj(_vfsproj),
        //_fsprojit(-2.0, 2.0), _vfsprojit(_fsprojit), 
        _chfsproj(_vfsproj),
        _calmetproj(_vfsproj), _chleproj(_vfsproj), _pvtxproj(),
        _svtxproj(_pvtxproj, _chfsproj, _jetaxes, 0.7, 2.0, 34e-3, 7.5, 34e-3)
    {

      setBeams(PROTON, ANTIPROTON);

      // Add particle/antiparticle vetoing: 12=nu_e, 14=nu_mu, 16=nu_tau
      _vfsproj
        .addVetoPairId(12)
        .addVetoPairId(14)
        .addVetoPairId(16);

      // Veto muons (PDG code = 13) with \f$ p_T \f$ above 1.0 GeV
      _vfsproj.addVetoDetail(13, 1.0, numeric_limits<double>::max());

      addProjection(_fsproj);
      addProjection(_vfsproj);
      addProjection(_chfsproj);
      addProjection(_jetsproj);
      addProjection(_calmetproj);
      addProjection(_chleproj);
      addProjection(_pvtxproj);
      addProjection(_svtxproj);
    }

    /// Factory method.
    static Analysis* create() { 
      return new CDF_2006_S6653332(); 
    }

    /// Get the name of this analysis.
    string getName() const {
      return "CDF_2006_S6653332";
    }

    /// Applying complex quality cuts on tracks to establish displaced vertices
    /// will be read as function pointer from the SVertex projection constructor.
    //static bool applyVtxTrackCuts(SVertex&, ParticleVector&, const HepMC::GenVertex& gpvtx, FourMomentum);

    void init();
    
    void analyze(const Event& event);
    
    void finalize();

  private:

    /// The final state projector.
    FinalState _fsproj;

    /// The visible final state projector.
    VetoedFinalState _vfsproj;

    /// The FastJets projector used by this analysis.
    FastJets _jetsproj;

    /// The charged final state projector.
    ChargedFinalState _chfsproj;

    /// The Calorimeter MET projector.
    TotalVisibleMomentum _calmetproj;

    /// The charged leptons projector.
    ChargedLeptons _chleproj;

    /// Projection to find the primary vertex.
    PVertex _pvtxproj;

    /// Projection to find secondary vertices.
    SVertex _svtxproj;
    

  private:

    /// Set of vectors defining the jet axes (for detached vertex finding).
    vector<FourMomentum> _jetaxes;
    
    /// Hide the assignment operator
    CDF_2006_S6653332& operator=(const CDF_2006_S6653332&);

    //@{
    /// Histograms
    AIDA::IHistogram1D* _histJetsPt;
    AIDA::IHistogram1D* _histJetsEta;
    AIDA::IHistogram1D* _histbJetsPt;
    AIDA::IHistogram1D* _histbJetsEta;
    //@}

  };

}

#endif
