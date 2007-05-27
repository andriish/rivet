// -*- C++ -*-
#ifndef RIVET_HepEx9506012_HH
#define RIVET_HepEx9506012_HH

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/FinalStateHCM.hh"
#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This analysis measures energy flow in DIS?
  class HepEx9506012 : public Analysis {

  public:

    /// The default constructor.
    inline HepEx9506012()
      : beams(), 
        lepton(beams, ELECTRON, POSITRON), 
        diskin(beams, lepton, PROTON),
        fsproj(lepton, diskin, fsp), 
        y1hcm(fsproj)
    { 
      setBeams(ELECTRON, PROTON);
      addProjection(beams);
      addProjection(lepton);
      addProjection(diskin);
      addProjection(fsproj);
      addProjection(y1hcm);
    }

  public:

    /// The name of this analysis is "HepEx9506012"
    inline string getName() const {
      return "HepEx9506012";
    }
    
    /// Initialize this analysis object.
    void init();
    
    /// Analyze one event.
    void analyze(const Event & event);
    
    /// Finalize this analysis object.
    void finalize();
    
    /// Return the RivetInfo object of this analysis object. 
    //RivetInfo getInfo() const;

  protected:

    /// Calculate the bin number from the DISKinematics projection.
    int getbin(const DISKinematics& dk);


  private:
    
    /// The Beam projector used.
    Beam beams;
    
    /// The DISLepton projector used.
    DISLepton lepton;
    
    /// The DISKinematics projector used.
    DISKinematics diskin;
    
    /// The FinalState projector used.
    FinalState fsp;
    
    /// The FinalStateHCM projector used.
    FinalStateHCM fsproj;
    
    /// The CentralEtHCM projector used.
    CentralEtHCM y1hcm;
    
    /// Some integer constants used.
    static const int nb = 24, nbin = 9;
    
    /// Some double constants used.
    static const double xmin, xmax;

    /// Histograms for the \f$ E_t \f$ flows
    vector<AIDA::IHistogram1D*> hEtFlow, hEtFlowStat;

    /// Histograms for averages in different kinematical bins.
    AIDA::IHistogram1D *hAvEt, *hAvX, *hAvQ2, *hN;

    /// Helper vector;
    vector<double> nev;

  private:

    /// Hidden assignment operator.
    HepEx9506012 & operator=(const HepEx9506012&);

  };

}

#endif
