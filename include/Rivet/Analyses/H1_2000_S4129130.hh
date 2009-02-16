// -*- C++ -*-
#ifndef RIVET_H1_2000_S4129130_HH
#define RIVET_H1_2000_S4129130_HH

#include "Rivet/Analysis.hh" 
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief H1 energy flow and charged particle spectra
  /// @author Peter Richardson
  /// Based on the HZtool analysis
  class H1_2000_S4129130 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    H1_2000_S4129130() {
      setBeams(ELECTRON, PROTON);
      const DISLepton& lepton = addProjection(DISLepton(ELECTRON, POSITRON), "Lepton");
      addProjection(DISKinematics(lepton, PROTON), "Kinematics");
      addProjection(FinalState(), "FS");
    }

    /// Factory method.
    static Analysis* create() { 
      return new H1_2000_S4129130(); 
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "4129130";
    }
    /// A short description of the analysis.
    string summary() const {
      return "H1 energy flow in DIS";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "H1";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2000";
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

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
    vector<AIDA::IHistogram1D *> _histETLowQa;
    vector<AIDA::IHistogram1D *> _histETHighQa;
    vector<AIDA::IHistogram1D *> _histETLowQb;
    vector<AIDA::IHistogram1D *> _histETHighQb;
    AIDA::IProfile1D * _histAverETCentral;
    AIDA::IProfile1D * _histAverETFrag;
    //@}

    /// @name storage of weights for normalisation
    //@{
    vector<double> _weightETLowQa;
    vector<double> _weightETHighQa;
    vector<double> _weightETLowQb;
    vector<double> _weightETHighQb;
    //@}
  };

}

#endif
