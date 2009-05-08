// -*- C++ -*-
#ifndef RIVET_ALEPH_1996_S3486095_HH
#define RIVET_ALEPH_1996_S3486095_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief ALEPH QCD study with event shapes and identified particles
  /// @author Holger Schulz
  class ALEPH_1996_S3486095 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ALEPH_1996_S3486095();

    /// Factory method.
    static Analysis* create() { 
      return new ALEPH_1996_S3486095(); 
    }

    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "3486095";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Studies of QCD with the ALEPH detector.";
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Summary paper of QCD results as measured by ALEPH at LEP 1. "
         << "The publication includes various event shape variables, "
         << "multiplicities (identified particles and inclusive), and "
         << "particle spectra.";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "ALEPH";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "LEP 1";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1996";
    }
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Holger Schulz <holger.schulz@physik.hu-berlin.de>";
      return ret;
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Hadronic Z decay events generated on the Z pole (sqrt(s) = 91.2 GeV)";
      return os.str();
    }

    string status() const {
      return "VALIDATED";
    }

    /// Journal, and preprint references
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys. Rept., 294, 1--165 (1998)";
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
    /// Little helper functions for the axis labels
    string unitdsigbyd(const string&);
    string unitdNbyd(const string&);
    string texmath(const string&);

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the 
    /// inclusive single particle distributions' normalisations.
    double _weightedTotalPartNum;
    double _weightedTotalNumPiPlus;       
    double _weightedTotalNumKPlus;      
    double _weightedTotalNumP;     
    double _weightedTotalNumPhoton;    
    double _weightedTotalNumPi0;   
    double _weightedTotalNumEta;  
    double _weightedTotalNumEtaPrime; 
    double _weightedTotalNumK0;
    double _weightedTotalNumLambda0;
    double _weightedTotalNumXiMinus;
    double _weightedTotalNumSigma1385Plus;
    double _weightedTotalNumXi1530_0;
    double _weightedTotalNumRho;
    double _weightedTotalNumOmega782;
    double _weightedTotalNumKStar892_0;
    double _weightedTotalNumPhi;
    double _weightedTotalNumKStar892Plus;
    double _numChParticles;

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_histSphericity;
    AIDA::IHistogram1D *_histAplanarity;

    AIDA::IHistogram1D *_hist1MinusT; 
    AIDA::IHistogram1D *_histTMinor; 
    
    AIDA::IHistogram1D *_histY3;
    AIDA::IHistogram1D *_histHeavyJetMass;
    AIDA::IHistogram1D *_histCParam;
    AIDA::IHistogram1D *_histOblateness; 
    
    AIDA::IHistogram1D *_histScaledMom; 
    AIDA::IHistogram1D *_histRapidityT;

    AIDA::IHistogram1D *_histPtSIn;
    AIDA::IHistogram1D *_histPtSOut;
    
    AIDA::IHistogram1D *_histJetRate2Durham;
    AIDA::IHistogram1D *_histJetRate3Durham;
    AIDA::IHistogram1D *_histJetRate4Durham;
    AIDA::IHistogram1D *_histJetRate5Durham;
   
    AIDA::IHistogram1D *_histLogScaledMom;
    
    
    AIDA::IHistogram1D *_histChMult;
    

    AIDA::IHistogram1D *_histMultiPiPlus;
    AIDA::IHistogram1D *_histMultiKPlus;
    AIDA::IHistogram1D *_histMultiP;
    AIDA::IHistogram1D *_histMultiPhoton;
    AIDA::IHistogram1D *_histMultiPi0;
    AIDA::IHistogram1D *_histMultiEta;
    AIDA::IHistogram1D *_histMultiEtaPrime;
    AIDA::IHistogram1D *_histMultiK0;
    AIDA::IHistogram1D *_histMultiLambda0;
    AIDA::IHistogram1D *_histMultiXiMinus;
    AIDA::IHistogram1D *_histMultiSigma1385Plus;
    AIDA::IHistogram1D *_histMultiXi1530_0;
    AIDA::IHistogram1D *_histMultiRho;
    AIDA::IHistogram1D *_histMultiOmega782;
    AIDA::IHistogram1D *_histMultiKStar892_0;
    AIDA::IHistogram1D *_histMultiPhi;
    AIDA::IHistogram1D *_histMultiKStar892Plus;
   
    // mean multiplicities
    AIDA::IHistogram1D *_histMeanChMult;
    AIDA::IHistogram1D *_histMeanChMultRapt05;
    AIDA::IHistogram1D *_histMeanChMultRapt10;
    AIDA::IHistogram1D *_histMeanChMultRapt15;
    AIDA::IHistogram1D *_histMeanChMultRapt20;
    
    AIDA::IHistogram1D *_histMeanMultiPi0;          
    AIDA::IHistogram1D *_histMeanMultiEta;          
    AIDA::IHistogram1D *_histMeanMultiEtaPrime;     
    AIDA::IHistogram1D *_histMeanMultiK0;           
    AIDA::IHistogram1D *_histMeanMultiRho;          
    AIDA::IHistogram1D *_histMeanMultiOmega782;        
    AIDA::IHistogram1D *_histMeanMultiPhi;         
    AIDA::IHistogram1D *_histMeanMultiKStar892Plus; 
    AIDA::IHistogram1D *_histMeanMultiKStar892_0;   
    AIDA::IHistogram1D *_histMeanMultiLambda0;      
    AIDA::IHistogram1D *_histMeanMultiSigma0;       
    AIDA::IHistogram1D *_histMeanMultiXiMinus;      
    AIDA::IHistogram1D *_histMeanMultiSigma1385Plus;
    AIDA::IHistogram1D *_histMeanMultiXi1530_0;     
    AIDA::IHistogram1D *_histMeanMultiOmegaOmegaBar;        
    //@}

  };

}

#endif
