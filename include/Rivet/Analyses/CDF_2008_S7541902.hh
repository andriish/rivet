// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7541902_HH
#define RIVET_CDF_2008_S7541902_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/SVertex.hh"

namespace Rivet {

  
  /// This CDF analysis provides jet pT distributions for 4 jet multplicity bins
  /// as well as the jet multiplicity distribution in W + jets events.
  /// e-Print: arXiv:0711.4044 [hep-ex]
  class CDF_2008_S7541902 : public Analysis {

  public:

    /// Constructor
    CDF_2008_S7541902() {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);

      // Veto neutrinos, W decay products, and muons with \f$ p_T \f$ above 1.0 GeV
      VetoedFinalState& vfs = *new VetoedFinalState();
      vfs
        .addVetoPairId(NU_E)
        .addVetoPairId(NU_MU)
        .addVetoPairId(NU_TAU)
        .addVetoDetail(13, 1.0, MAXDOUBLE) //< @todo Check that this also covers antimuons
        .addDecayProductsVeto(WPLUSBOSON)
        .addDecayProductsVeto(WMINUSBOSON);
      addProjection(vfs, "FS");
      addProjection(*new FastJets(getProjection<FinalState>("FS"),
                                  FastJets::CDFJETCLU, 0.4), "Jets");
    }


    /// Factory method.
    static Analysis* create() { 
      return new CDF_2008_S7541902(); 
    }


    /// @name Publication metadata
    //@{
    /// Get the SPIRES ID
    string getSpiresId() const {
      return "7541902";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Jet pT distributions for 4 jet multiplicity bins as well as the jet multiplicity distribution in W + jets events.";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2008";
    }
    /// Journal, and preprint references.
    vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("arXiv:0711.4044 [hep-ex]");
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
    
    //@{
    /// Histograms
    AIDA::IHistogram1D* _histJetMult;
    AIDA::IHistogram1D* _histJetEt[4];
    //@}

    /// Hide the assignment operator
    CDF_2008_S7541902& operator=(const CDF_2008_S7541902&);
  };

}

#endif
