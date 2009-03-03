// -*- C++ -*-
#ifndef RIVET_CDF_2008_S8095620_HH
#define RIVET_CDF_2008_S8095620_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
namespace Rivet {


  /// Implementation of CDF Run II Z+b-jet cross section paper

  class CDF_2008_S8095620 : public Analysis {

  public:

    /// Constructor.

    /// jet cuts: |eta| <= 1.5


    CDF_2008_S8095620()
      : _Rjet(0.7), _JetPtCut(20.), _JetEtaCut(1.5)
    { 
      setBeams(PROTON, ANTIPROTON);

      const FinalState fs(-3.6, 3.6);
      addProjection(fs, "FS");

      // Create a final state with any e+e- or mu+mu- pair with 
      // invariant mass 76 -> 106 GeV and ET > 20 (Z decay products)
      std::vector<std::pair<long,long> > vids;
      vids.push_back(make_pair(ELECTRON, POSITRON));
      vids.push_back(make_pair(MUON, ANTIMUON));
      FinalState fs2(-3.6, 3.6, 20*GeV);
      InvMassFinalState invfs(fs2, vids, 76*GeV, 106*GeV);
      addProjection(invfs, "INVFS");
      // Make a final state without the Z decay products for jet clustering
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(invfs);
      addProjection(vfs, "VFS");
      addProjection(FastJets(vfs, FastJets::CDFMIDPOINT, 0.7), "Jets");

    }


    /// Factory method
    static Analysis* create() { 
      return new CDF_2008_S8095620(); 
    }


    /// @name Publication metadata
    //@{
    /// Get SPIRES ID code.
    string spiresId() const {
      return "7782535";
    }
    /// A short description of the analysis.
    string summary() const {
      return "CDF Run II Z+b-jet cross section paper, 2 fb-1";
    }

/// A full description of the analysis.
    string description() const {
      ostringstream os;
      os <<  "Measurement of the b-jet production cross section for " 
	 <<  "events containing a $Z$ boson produced in $p\\bar{p}$ collisions at "
	 <<  "$\\sqrt{s}=1.96$ TeV, using data corresponding to an integrated "
	 <<  "luminosity of 2 fb$^{-1}$ collected by the CDF II detector at the "
	 <<  "Tevatron. $Z$ bosons are selected in the electron and muon decay "
	 <<  "modes. Jets are considered with transverse energy $E_T>20$ GeV and "
	 <<  "pseudorapidity $|\\eta|<1.5$. The ratio of the integrated $Z$ + b-jet "
	 <<  "cross section to the inclusive $Z$ production cross section is "
	 <<  "measured differentially in jet $E_T$, jet $\\eta$, "
	 <<  "$Z$-boson transverse momentum, number of jets, and number of b-jets. "
	 <<  "The first two measurements have an entry for each b-jet in the event, "
	 <<  "the last three measurements have one entry per event.";
      return os.str();
    }
  /// Characteristics of events to be processed by this analysis
    string runInfo() const {
      ostringstream os;
      os << "Requires the process $p\\bar{p} \\rightarrow {Z} "
         <<  "\\rightarrow{\\ell}\\ell$, where $\\ell$ is $e$ or $\\mu$. "
         << "Additional hard jets will also have to be included to get a good description. ";
      return os.str();
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }

  /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 2";
    }

    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2008";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Emily Nurse <nurse@hep.ucl.ac.uk>";
      return rtn;
    }

    /// Journal, and preprint references.
    virtual vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0812.4458");
      return ret;
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}


  private:

    double _Rjet;
    double _JetPtCut;
    double _JetEtaCut;
    //@{
    /// Histograms
    AIDA::IHistogram1D* _dSdET;
    AIDA::IHistogram1D* _dSdETA;
    AIDA::IHistogram1D* _dSdNJet; 
    AIDA::IHistogram1D* _dSdNbJet; 
    AIDA::IHistogram1D* _dSdZpT; 

    //@}

  };

}

#endif
