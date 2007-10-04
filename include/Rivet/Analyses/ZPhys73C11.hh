// -*- C++ -*-
#ifndef RIVET_ZPHYS73C11_HH
#define RIVET_ZPHYS73C11_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {


  /// Implementation of DELPHI event shape paper
  class ZPhys73C11 : public Analysis {

  public:

    /// Default constructor.
    inline ZPhys73C11()
      : _cnfsproj(), _cfsproj(_cnfsproj),
        _cspherproj(_cfsproj), _cnspherproj(_cnfsproj), 
        _thrustproj(_cfsproj), 
        _parisiproj(_cfsproj),
        _hemiproj(_cfsproj)
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(_cfsproj);
      addProjection(_cnfsproj);
      addProjection(_beamsproj);
      addProjection(_cspherproj);
      addProjection(_cnspherproj);
      addProjection(_thrustproj);
      addProjection(_parisiproj);
      addProjection(_hemiproj);
    }

    /// Factory method.
    static Analysis* create() { 
      return new ZPhys73C11(); 
    }

    /// Get the name of this analysis.
    inline string getName() const {
      return "ZPhys73C11";
    }

    virtual void init();

    virtual void analyze(const Event& event);

    virtual void finalize();


  private:

    /// Hide the assignment operator
    ZPhys73C11& operator=(const ZPhys73C11&);


  private:

    /// The final state projector.
    FinalState _cnfsproj;

    /// Charged final state projector.
    ChargedFinalState _cfsproj;

    /// Projection to get the beams.
    Beam _beamsproj;

    /// Sphericity projection for charged particles.
    Sphericity _cspherproj;

    /// Sphericity projection for all particles.
    Sphericity _cnspherproj;

    /// Thrust projection.
    Thrust _thrustproj;

    /// Parisi tensor (C and D params) projection.
    ParisiTensor _parisiproj;

    /// Projection to calculate event hemisphere masses and broadenings.
    Hemispheres _hemiproj;


    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_histPtTIn,  *_histPtTInCN;
    AIDA::IHistogram1D *_histPtTOut, *_histPtTOutCN;
    AIDA::IHistogram1D *_histPtSIn,  *_histPtSInCN;
    AIDA::IHistogram1D *_histPtSOut, *_histPtSOutCN;

    AIDA::IHistogram1D *_histRapidityT, *_histRapidityTCN;
    AIDA::IHistogram1D *_histRapidityS, *_histRapiditySCN;

    AIDA::IHistogram1D *_histScaledMom, *_histLogScaledMom;

    AIDA::IHistogram1D *_histPtTOutVsXp, *_histPtTInVsXp;

    AIDA::IHistogram1D *_hist1MinusT, *_hist1MinusTCN; 
    AIDA::IHistogram1D *_histTMajor, *_histTMajorCN; 
    AIDA::IHistogram1D *_histTMinor, *_histTMinorCN; 
    AIDA::IHistogram1D *_histOblateness, *_histOblatenessCN; 

    AIDA::IHistogram1D *_histSphericity, *_histSphericityCN;
    AIDA::IHistogram1D *_histAplanarity, *_histAplanarityCN;
    AIDA::IHistogram1D *_histPlanarity, *_histPlanarityCN;

    AIDA::IHistogram1D *_histCParam, *_histCParamCN;
    AIDA::IHistogram1D *_histDParam, *_histDParamCN;

    AIDA::IHistogram1D *_histHemiMassD, *_histHemiMassDCN;
    AIDA::IHistogram1D *_histHemiMassH, *_histHemiMassHCN;
    AIDA::IHistogram1D *_histHemiMassL, *_histHemiMassLCN;
               
    AIDA::IHistogram1D *_histHemiBroadW, *_histHemiBroadWCN;
    AIDA::IHistogram1D *_histHemiBroadN, *_histHemiBroadNCN;
    AIDA::IHistogram1D *_histHemiBroadT, *_histHemiBroadTCN;
    AIDA::IHistogram1D *_histHemiBroadD, *_histHemiBroadDCN;

    AIDA::IHistogram1D *_histDiffRate2Durham, *_histDiffRate2DurhamCN;
    AIDA::IHistogram1D *_histDiffRate2Jade,   *_histDiffRate2JadeCN; 
    AIDA::IHistogram1D *_histDiffRate3Durham, *_histDiffRate3DurhamCN;
    AIDA::IHistogram1D *_histDiffRate3Jade,   *_histDiffRate3JadeCN;
    AIDA::IHistogram1D *_histDiffRate4Durham, *_histDiffRate4DurhamCN;
    AIDA::IHistogram1D *_histDiffRate4Jade,   *_histDiffRate4JadeCN;

    AIDA::IHistogram1D *_histEEC, *_histAEEC;
    //@}
  
  };

}

#endif
