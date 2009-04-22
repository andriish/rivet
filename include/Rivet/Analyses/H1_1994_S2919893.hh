// -*- C++ -*-
#ifndef RIVET_H1_1994_S2919893_HH
#define RIVET_H1_1994_S2919893_HH

#include "Rivet/Analysis.hh" 

namespace Rivet {


  /// @brief H1 energy flow and charged particle spectra
  /// @author Peter Richardson
  /// Based on the equivalent HZTool analysis
  class H1_1994_S2919893 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    H1_1994_S2919893();

    /// Factory method.
    static Analysis* create() { 
      return new H1_1994_S2919893(); 
    }

    //@}
    

    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "2919893";
    }
    /// A short description of the analysis.
    string summary() const {
      return "H1 energy flow and charged particle spectra in DIS";
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Global properties of the hadronic final state in deep inelastic scattering events "
         << "at HERA are investigated. The data are corrected for detector effects. Energy flows "
         << "in both the laboratory frame and the hadronic centre of mass system, and "
         << "energy-energy correlations in the laboratory frame are presented."
         << "\n\n"
         << "Historically, the Ariadne colour dipole model provided the only satisfactory description "
         << "of this data, hence making it a useful 'target' analysis for MC shower models.";
      return os.str();
    }
    /// Event type required by this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Event type: e- p / e+ p (TODO: which?) deep inelastic scattering\n"
         << "* Energy: ???";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "H1";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "HERA";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1994";
    }
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Peter Richardson <peter.richardson@durham.ac.uk>";
      return rtn;
    }
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Z.Phys.C63:377-390,1994";
      ret += "doi:10.1007/BF01580319";
      return ret;
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
