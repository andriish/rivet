// -*- C++ -*-
#ifndef RIVET_CDF_2006_S6653332_HH
#define RIVET_CDF_2006_S6653332_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {


  /* @brief CDF Run II analysis: jet \f$ p_T \f$ and \f$ \eta \f$ distributions in Z + (b) jet production
   * @author Lars Sonnenschein
   *
   * This CDF analysis provides \f$ p_T \f$ and \f$ \eta \f$ distributions of
   * jets in Z + (b) jet production, before and after tagging.
   */
  class CDF_2006_S6653332 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{
    /// Constructor
    CDF_2006_S6653332()  : _Rjet(0.7), _JetPtCut(20.), _JetEtaCut(1.5),
			   _sumWeightsWithZ(0.0), _sumWeightsWithZJet(0.0)
    { 
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
      const FinalState fs(-3.6, 3.6);
      addProjection(fs, "FS");

      // Create a final state with any e+e- or mu+mu- pair with 
      // invariant mass 76 -> 106 GeV and ET > 20 (Z decay products)
      std::vector<std::pair<long,long> > vids;
      vids.push_back(make_pair(ELECTRON, POSITRON));
      vids.push_back(make_pair(MUON, ANTIMUON));
      FinalState fs2(-3.6, 3.6);
      InvMassFinalState invfs(fs2, vids, 76*GeV, 106*GeV);
      addProjection(invfs, "INVFS");
      // Make a final state without the Z decay products for jet clustering
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(invfs);
      addProjection(vfs, "VFS");
      addProjection(FastJets(vfs, FastJets::CDFMIDPOINT, 0.7), "Jets");

    }


    /// Factory method.
    static Analysis* create() { 
      return new CDF_2006_S6653332(); 
    }
    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "6653332";
    }

    /// A short description of the analysis.
    string summary() const {
      return "pT and eta distributions of jets in Z + jet production";
    }

    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Measurement of the b jet cross section in events with Z boson in p pbar "
         << "collisions at center-of-mass energy sqrt(s) = 1.96 TeV. The data cover jet "
         << "transverse momenta above 20 GeV and jet pseudo-rapidities in the range "
         << "-1.5 to 1.5. Z bosons are identified in their electron and muon decay modes "
         << "in an invariant dilepton mass range between 66 and 116 GeV."
        //   << "\n\n"
         << "";
      return os.str();
    }

    /// Event type required by this analysis.
    string runInfo() const {
      ostringstream os;
      os << "* Energy: sqrt(s) = 1960 GeV\n"
         << "* Event type: Z + jets events\n"
         << "* Jets min pT cut: pT_{jet} > 20 GeV\n"
         << "* Leptons min pT cut: pT_{jet} > 10 GeV";
      return os.str();
    }

    /// Experiment which performed and published this analysis. 
    string experiment() const {
      return "CDF";
    }

    /// Collider on which the experiment ran
    string collider() const {
      return "Tevatron Run 2";
    }

    /// When published (preprint year according to SPIRES). 
    string year() const {
      return "2006";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Lars Sonnenschein <Lars.Sonnenschein@cern.ch>";
      return rtn;
    }

    /// Publication references.
    vector<string> references() const {
      vector<string> ret;
      ret += "Phys.Rev.D.74:032008,2006";
      ret += "doi:10.1103/PhysRevD.74.032008";
      ret += "arXiv:hep-ex/0605099v2";
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


    double _Rjet;
    double _JetPtCut;
    double _JetEtaCut;

    double _sumWeightsWithZ;
    double _sumWeightsWithZJet;


    //@{
    /// Histograms
    AIDA::IHistogram1D* _sigmaBJet;
    AIDA::IHistogram1D* _ratioBJetToZ;
    AIDA::IHistogram1D* _ratioBJetToJet;
    //@}

  };


}

#endif
