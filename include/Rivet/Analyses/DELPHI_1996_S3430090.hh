// -*- C++ -*-
#ifndef RIVET_DELPHI_1996_S3430090_HH
#define RIVET_DELPHI_1996_S3430090_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {


  /// Implementation of DELPHI event shape paper
  class DELPHI_1996_S3430090 : public Analysis {

  public:

    /// Default constructor.
    DELPHI_1996_S3430090()
      : _cnfsproj(), _cfsproj(_cnfsproj),
#ifdef __HAVE_JADE
        _cjadejetproj(_cfsproj, FastJets::JADE, 0.7),
        _cnjadejetproj(_cnfsproj, FastJets::JADE, 0.7),
        _cdurjetproj(_cfsproj, FastJets::DURHAM, 0.7),
        _cndurjetproj(_cnfsproj, FastJets::DURHAM, 0.7),
#else
        _cjadejetproj(_cfsproj, FastJets::KT, 0.7),
        _cnjadejetproj(_cnfsproj, FastJets::KT, 0.7),
        _cdurjetproj(_cfsproj, FastJets::KT, 0.7),
        _cndurjetproj(_cnfsproj, FastJets::KT, 0.7),
#endif
        _cspherproj(_cfsproj), _cnspherproj(_cnfsproj), 
        _cthrustproj(_cfsproj), _cnthrustproj(_cnfsproj), 
        _cparisiproj(_cfsproj), _cnparisiproj(_cnfsproj),
        _chemiproj(_cfsproj, _cthrustproj), 
        _cnhemiproj(_cnfsproj, _cnthrustproj)
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(_cfsproj);
      addProjection(_cnfsproj);
      addProjection(_beamsproj);
      addProjection(_cspherproj);
      addProjection(_cnspherproj);
      addProjection(_cthrustproj);
      addProjection(_cnthrustproj);
      addProjection(_cparisiproj);
      addProjection(_cnparisiproj);
      addProjection(_chemiproj);
      addProjection(_cnhemiproj);
    }


    /// Factory method.
    static Analysis* create() { 
      return new DELPHI_1996_S3430090(); 
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "3430090";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Delphi MC tuning on event shapes and identified particles.";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "DELPHI";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "1996";
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /// Hide the assignment operator
    DELPHI_1996_S3430090& operator=(const DELPHI_1996_S3430090&);

    /// The final state projector.
    FinalState _cnfsproj;

    /// Charged final state projector.
    ChargedFinalState _cfsproj;

    /// Projection to get the beams.
    Beam _beamsproj;

    /// Jet algorithms
    FastJets _cjadejetproj, _cnjadejetproj;
    FastJets _cdurjetproj, _cndurjetproj;

    /// Sphericity projections.
    Sphericity _cspherproj, _cnspherproj;

    /// Thrust projections.
    Thrust _cthrustproj, _cnthrustproj;

    /// Parisi tensor (C and D params) projections.
    ParisiTensor _cparisiproj, _cnparisiproj;

    /// Projections to calculate event hemisphere masses and broadenings.
    Hemispheres _chemiproj, _cnhemiproj;

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the 
    /// inclusive single particle distributions' normalisations.
    double _weightedTotalPartNumC, _weightedTotalPartNumCN;
    
    /// Store the sum of weights for events which are rejected by the
    /// lepton veto. We need this number to get the normalisation of single
    /// particle distributions right.
    double _sumOfRejectedWeights;

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_histPtTInC,  *_histPtTInCN;
    AIDA::IHistogram1D *_histPtTOutC, *_histPtTOutCN;
    AIDA::IHistogram1D *_histPtSInC,  *_histPtSInCN;
    AIDA::IHistogram1D *_histPtSOutC, *_histPtSOutCN;

    AIDA::IHistogram1D *_histRapidityTC, *_histRapidityTCN;
    AIDA::IHistogram1D *_histRapiditySC, *_histRapiditySCN;

    AIDA::IHistogram1D *_histScaledMom, *_histLogScaledMom;

    AIDA::IProfile1D   *_histPtTOutVsXp, *_histPtVsXp;

    AIDA::IHistogram1D *_hist1MinusTC, *_hist1MinusTCN; 
    AIDA::IHistogram1D *_histTMajorC, *_histTMajorCN; 
    AIDA::IHistogram1D *_histTMinorC, *_histTMinorCN; 
    AIDA::IHistogram1D *_histOblatenessC, *_histOblatenessCN; 

    AIDA::IHistogram1D *_histSphericityC, *_histSphericityCN;
    AIDA::IHistogram1D *_histAplanarityC, *_histAplanarityCN;
    AIDA::IHistogram1D *_histPlanarityC, *_histPlanarityCN;

    AIDA::IHistogram1D *_histCParamC, *_histCParamCN;
    AIDA::IHistogram1D *_histDParamC, *_histDParamCN;

    AIDA::IHistogram1D *_histHemiMassDC, *_histHemiMassDCN;
    AIDA::IHistogram1D *_histHemiMassHC, *_histHemiMassHCN;
    AIDA::IHistogram1D *_histHemiMassLC, *_histHemiMassLCN;
               
    AIDA::IHistogram1D *_histHemiBroadWC, *_histHemiBroadWCN;
    AIDA::IHistogram1D *_histHemiBroadNC, *_histHemiBroadNCN;
    AIDA::IHistogram1D *_histHemiBroadTC, *_histHemiBroadTCN;
    AIDA::IHistogram1D *_histHemiBroadDC, *_histHemiBroadDCN;

    AIDA::IHistogram1D *_histDiffRate2DurhamC, *_histDiffRate2DurhamCN;
    AIDA::IHistogram1D *_histDiffRate2JadeC,   *_histDiffRate2JadeCN; 
    AIDA::IHistogram1D *_histDiffRate3DurhamC, *_histDiffRate3DurhamCN;
    AIDA::IHistogram1D *_histDiffRate3JadeC,   *_histDiffRate3JadeCN;
    AIDA::IHistogram1D *_histDiffRate4DurhamC, *_histDiffRate4DurhamCN;
    AIDA::IHistogram1D *_histDiffRate4JadeC,   *_histDiffRate4JadeCN;

    AIDA::IHistogram1D *_histEEC, *_histAEEC;
    //@}
  
  };

}

#endif
