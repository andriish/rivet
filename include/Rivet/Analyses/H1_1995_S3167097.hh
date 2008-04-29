// -*- C++ -*-
#ifndef RIVET_H1_1995_S3167097_HH
#define RIVET_H1_1995_S3167097_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalStateHCM.hh"
#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This analysis measures energy flow in DIS?
  /// @author Leif Lonnblad
  class H1_1995_S3167097 : public Analysis {
  public:

    /// Constructor.
    H1_1995_S3167097() { 
      setBeams(ELECTRON, PROTON);
      const DISLepton& lepton = addProjection(*new DISLepton(ELECTRON, POSITRON), "Lepton");
      const DISKinematics& diskin = addProjection(*new DISKinematics(lepton, PROTON), "Kinematics");
      const FinalStateHCM& fshcm = addProjection(*new FinalStateHCM(diskin), "FS");
      addProjection(*new CentralEtHCM(fshcm), "Y1HCM");
      addCut("x", MORE_EQ, _xmin);
      addCut("x", LESS_EQ, _xmax);
    }

  public:

    /// Factory method
    static Analysis* create() { return new H1_1995_S3167097(); }

    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "3167097";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "TODO";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "H1";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "1995";
    }
    //@}
    

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

    
  private:

    /// Calculate the bin number from the DISKinematics projection.
    int _getbin(const DISKinematics& dk);


    /// Some integer constants used.
    static const size_t _nb = 24, _nbin = 9;
    
    /// Some double constants used.
    static const double _xmin, _xmax;

    /// Histograms for the \f$ E_T \f$ flows
    vector<AIDA::IHistogram1D*> _hEtFlow, _hEtFlowStat;

    /// Histograms for averages in different kinematical bins.
    AIDA::IHistogram1D *_hAvEt, *_hAvX, *_hAvQ2, *_hN;

    /// Helper vector;
    vector<double> _nev;

    /// Hidden assignment operator.
    H1_1995_S3167097& operator=(const H1_1995_S3167097&);

  };

}

#endif
