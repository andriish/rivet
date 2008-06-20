// -*- C++ -*-
#ifndef RIVET_H1_1994_S2919893_HH
#define RIVET_H1_1994_S2919893_HH

#include "Rivet/Analysis.hh" 
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// Implementation of H1 energy flow and charged particle spectra
  /// @author Peter Richardson
  /// based on the equivalent HZTool analysis
  class H1_1994_S2919893 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    H1_1994_S2919893() {
      setBeams(ELECTRON, PROTON);
      const DISLepton& lepton = 
	addProjection(DISLepton(ELECTRON, POSITRON), "Lepton");
      addProjection(DISKinematics(lepton, PROTON), "Kinematics");
      addProjection(FinalState(), "FS");
    }


    /// Factory method.
    static Analysis* create() { 
      return new H1_1994_S2919893(); 
    }

    //@}


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "2919893";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "H1 energy flow and charged particle spectra in DIS";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "H1";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "1994";
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
    H1_1994_S2919893& operator=(const H1_1994_S2919893&);

    /**
     *  Polar angle with right direction of the beam
     */
    inline double beamAngle(const FourVector& v, const bool & order) {
      double thel = v.polarAngle()/degree;
      if(thel<0.) thel+=180.;
      if(!order) thel = 180.-thel;
      return thel;
    }

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_histEnergyFlowLowX;
    AIDA::IHistogram1D *_histEnergyFlowHighX;
    AIDA::IHistogram1D *_histEECLowX;
    AIDA::IHistogram1D *_histEECHighX;
    AIDA::IHistogram1D *_histSpectraW77;
    AIDA::IHistogram1D *_histSpectraW122;
    AIDA::IHistogram1D *_histSpectraW169;
    AIDA::IHistogram1D *_histSpectraW117;
    AIDA::IProfile1D *_histPT2;
    //@}

    /// @name storage of weight to calculate averages for normalisation
    //@{
    pair<double,double> _w77,_w122,_w169,_w117,_wEnergy;
    //@}
  };

}

#endif
