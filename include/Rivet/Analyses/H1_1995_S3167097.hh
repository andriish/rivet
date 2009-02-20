// -*- C++ -*-
#ifndef RIVET_H1_1995_S3167097_HH
#define RIVET_H1_1995_S3167097_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalStateHCM.hh"
#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Measures energy flow in DIS? To be checked!
  /// @todo Check this analysis!
  /// @author Leif Lonnblad
  class H1_1995_S3167097 : public Analysis {
  public:

    /// Constructor.
    H1_1995_S3167097() { 
      setBeams(ELECTRON, PROTON);
      const DISLepton& lepton = addProjection(DISLepton(ELECTRON, POSITRON), "Lepton");
      const DISKinematics& diskin = addProjection(DISKinematics(lepton, PROTON), "Kinematics");
      const FinalStateHCM& fshcm = addProjection(FinalStateHCM(diskin), "FS");
      addProjection(CentralEtHCM(fshcm), "Y1HCM");
      //addCut("x", MORE_EQ, _xmin);
      //addCut("x", LESS_EQ, _xmax);
    }

  public:

    /// Factory method
    static Analysis* create() { return new H1_1995_S3167097(); }

    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "3167097";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Transverse energy and forward jet production in the low x regime at H1";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "H1";
    }
    
    /// Collider on which the experiment ran.
    string collider() const {
      return "HERA Run I";
    }
    
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1995";
    }
    
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn.push_back("Leif Lonnblad <leif.lonnblad@thep.lu.se>");
      return rtn;
    }
    
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "DIS events at low x may be sensitive to new QCD dynamics such as BFKL or CCFM radiation."
         << " In particular, BFKL is expected to produce more radiation at high transverse energy "
         << " in the rapidity span between the proton remnant and the struck quark jet." 
         << " Performing a transverse energy sum in bins of x and $\\eta$ may distinguish" 
         << " between DGLAP and BFKL evolution.";
      return os.str();
    }
    
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "HERA beam conditions: 820~GeV protons colliding with 26.7~GeV electrons."
         << "DIS events with an outgoing electron energy > 12~GeV."
         << "5~GeV$^2 < Q^2 <$100~GeV$^2$, $10^{-4} < x < 10^{-2}$.";
      return os.str();
    }
    
    /// Status of this routine (VALIDATED or UNVALIDATED)
    string status() const {
      return "UNVALIDATED";
    }
    
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("Phys.Lett.B356:118,1995");
      ret.push_back("hep-ex/9506012");
      return ret;
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

  };

}

#endif
