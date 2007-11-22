// -*- C++ -*-
#ifndef RIVET_S3167097_HH
#define RIVET_S3167097_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalStateHCM.hh"
#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This analysis measures energy flow in DIS?
  class S3167097 : public Analysis {

  public:

    /// The default constructor.
    inline S3167097()
      : _fsproj(), _beamsproj(), 
        _leptonproj(_beamsproj, _fsproj, ELECTRON, POSITRON), 
        _diskinproj(_beamsproj, _leptonproj, PROTON),
        _fshcmproj(_leptonproj, _diskinproj, _fsproj), 
        _y1hcmproj(_fshcmproj)
    { 
      setBeams(ELECTRON, PROTON);
      addProjection(_beamsproj);
      addProjection(_leptonproj);
      addProjection(_diskinproj);
      addProjection(_fshcmproj);
      addProjection(_y1hcmproj);
      addCut("x", MORE_EQ, _xmin);
      addCut("x", LESS_EQ, _xmax);
    }

  public:

    /// Factory method
    static Analysis* create() { return new S3167097(); }

    /// Get the name of this analysis.
    inline string getName() const {
      return "S3167097";
    }
    
    /// Initialize this analysis object.
    void init();
    
    /// Analyze one event.
    void analyze(const Event& event);
    
    /// Finalize this analysis object.
    void finalize();
    
  protected:

    /// Calculate the bin number from the DISKinematics projection.
    int getbin(const DISKinematics& dk);


  private:

    /// The FinalState projector used.
    FinalState _fsproj;
    
    /// The Beam projector used.
    Beam _beamsproj;
    
    /// The DISLepton projector used.
    DISLepton _leptonproj;
    
    /// The DISKinematics projector used.
    DISKinematics _diskinproj;
    
    /// The FinalStateHCM projector used.
    FinalStateHCM _fshcmproj;
    
    /// The CentralEtHCM projector used.
    CentralEtHCM _y1hcmproj;
    
    /// Some integer constants used.
    static const int _nb = 24, _nbin = 9;
    
    /// Some double constants used.
    static const double _xmin, _xmax;

    /// Histograms for the \f$ E_T \f$ flows
    vector<AIDA::IHistogram1D*> _hEtFlow, _hEtFlowStat;

    /// Histograms for averages in different kinematical bins.
    AIDA::IHistogram1D *_hAvEt, *_hAvX, *_hAvQ2, *_hN;

    /// Helper vector;
    vector<double> _nev;

  private:

    /// Hidden assignment operator.
    S3167097& operator=(const S3167097&);

  };

}

#endif
