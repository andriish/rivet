// -*- C++ -*-
#ifndef RIVET_ZPHYS73C11_HH
#define RIVET_ZPHYS73C11_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {


  /// Implementation of DELPHI event shape paper
  class ZPhys73C11 : public Analysis {

  public:

    /// Default constructor.
    inline ZPhys73C11()
      : _spherproj(_fsproj)
    { 
      setBeams(ELECTRON, POSITRON); 
      addProjection(_fsproj);
      addProjection(_beamsproj);
      addProjection(_spherproj);
    }

    /// Factory method
    static Analysis* create() { 
      return new ZPhys73C11(); 
    }

    /// The name of this analysis is "ZPhys73C11"
    inline string getName() const {
      return "ZPhys73C11";
    }

    virtual void init();

    virtual void analyze(const Event& event);

    virtual void finalize();

  private:

    /// The projectors used by this analysis.
    FinalState _fsproj;

    Beam _beamsproj;

    Sphericity _spherproj;

  private:

    /// Hide the assignment operator
    ZPhys73C11& operator=(const ZPhys73C11&);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histLogScaledMom;
    AIDA::IHistogram1D* _histScaledMom;
    AIDA::IHistogram1D* _histRapidityT;
    AIDA::IHistogram1D* _histRapidityS;
    AIDA::IHistogram1D* _histPtTIn;
    AIDA::IHistogram1D* _histPtTOut;
    AIDA::IHistogram1D* _histPtTInVsXp;
    AIDA::IHistogram1D* _histPtTOutVsXp;

    AIDA::IHistogram1D* _hist1MinusT; 
    AIDA::IHistogram1D* _histTMajor; 
    AIDA::IHistogram1D* _histTMinor; 
    AIDA::IHistogram1D* _histOblateness; 

    AIDA::IHistogram1D* _histDiffRate2Durham;
    AIDA::IHistogram1D* _histDiffRate2Jade; 
    AIDA::IHistogram1D* _histDiffRate2JadeCN;
    AIDA::IHistogram1D* _histDiffRate3Durham;
    AIDA::IHistogram1D* _histDiffRate3Jade;
    AIDA::IHistogram1D* _histDiffRate3JadeCN;
    AIDA::IHistogram1D* _histDiffRate4Durham;
    AIDA::IHistogram1D* _histDiffRate4Jade;

    AIDA::IHistogram1D* _histSphericity;
    AIDA::IHistogram1D* _histAplanarity;
    AIDA::IHistogram1D* _histPlanarity;
    AIDA::IHistogram1D* _histSphericityCN; 
    AIDA::IHistogram1D* _histAplanarityCN;

    AIDA::IHistogram1D* _histHemiMassD;
    AIDA::IHistogram1D* _histHemiMassH;
    AIDA::IHistogram1D* _histHemiMassL;
               
    AIDA::IHistogram1D* _histHemiBroadW;
    AIDA::IHistogram1D* _histHemiBroadN;
    AIDA::IHistogram1D* _histHemiBroadT;
    AIDA::IHistogram1D* _histHemiBroadD;

    AIDA::IHistogram1D* _histCParam;
    AIDA::IHistogram1D* _histDParam;
           
    AIDA::IHistogram1D* _histEEC;
    AIDA::IHistogram1D* _histAEEC;
  //@}
  
  };

}

#endif
