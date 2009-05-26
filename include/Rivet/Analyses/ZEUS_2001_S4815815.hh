// -*- C++ -*-
#ifndef RIVET_ZEUS_2001_S4815815_HH
#define RIVET_ZEUS_2001_S4815815_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

  /// @brief ZEUS dijet photoproduction study used in the ZEUS Jets PDF fit
  ///
  /// This class is a reproduction of the HZTool routine for the ZEUS 
  /// dijet photoproduction paper which was used in the ZEUS Jets PDF fit.  
  ///
  /// @author Jon Butterworth
  class ZEUS_2001_S4815815 : public Analysis {

  public:

    /// Default constructor.
    ZEUS_2001_S4815815();

    /// Factory method.
    static Analysis* create() { 
      return new ZEUS_2001_S4815815(); 
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string spiresId() const {
      return "4815815";
    }
    /// Get a description of the analysis.
    string summary() const {
      return "Dijet photoproduction analysis";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "ZEUS";
    }
    /// Collider on which the experiment ran.
    string collider() const{
      return "HERA Run I";
    }
    
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2001";
    }
    
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn.push_back("Jon Butterworth <jmb@hep.ucl.ac.uk>");
      return rtn;
    }
    
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os<< "ZEUS photoproduction of jets from proton-positron collisions at beam energies of 820~GeV on 27.5~GeV."
        << " Photoproduction can either be direct, in which case the photon interacts directly with the parton," 
        << " or resolved, in which case the photon acts as a source of quarks and gluons."
        << " A photon-proton centre of mass energy of between 134~GeV and 227~GeV is probed," 
        << " with values of xP, the fractional momentum of the partons inside the proton, predominantly in the region"
        << " between 0.01 and 0.1. The fractional momentum of the partons from the photon, $x\\gamma$, is in the region 0.1 to 1."
        << " Jets are reconstructed in the range $-1<|\\eta|<2.4$ using the kT algorithm with an R parameter of 1.0."
        << " The minimum pT of the leading jet should be greater then 14~GeV, and at least one other jet must have pT>11~GeV.";
      return os.str();
    }
    
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* HERA beam conditions: 820~GeV protons colliding with 27.5~GeV positrons\n"
         << "* Direct and resolved photoproduction of di-jets\n"
         << "* Leading jet pT > 14~GeV, second jet pT > 11~GeV\n"
         << "* Jet pseudorapidity $-1 < \\eta < 2.4$";
      return os.str();
    }
    
    /// Status of this routine (VALIDATED or UNVALIDATED)
    string status() const{
      return "UNVALIDATED";
    }
    
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> refs;
      refs.push_back( "Eur.Phys.J.C23:615,2002" );
      refs.push_back( "DESY 01/220" );
      refs.push_back("hep-ex/0112029");
      return refs;
    }
    
    //@}

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histJetEt1;
    //@}

  };


}

#endif
