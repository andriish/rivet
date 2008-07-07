// -*- C++ -*-
#ifndef RIVET_CDF_2005_S6217184_HH
#define RIVET_CDF_2005_S6217184_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// @brief CDF Run II jet shape paper 
  /// e-print: arXiv:hep-ex/0505013 
  class CDF_2005_S6217184 : public Analysis {

  public:

    /// Constructor.
    /// \f$ \eta \in [-2,2] \f$ cut used on final state.
    /// Jet shape \f$ r_\text{min} = 0.0 \f$, \f$ r_\text{max} = 0.7 \f$, 
    /// interval = 0.1, r1minPsi = 0.3.
    CDF_2005_S6217184()
      : _pvzmax(600*mm), _Rjet(0.7)
    { 
      setBeams(PROTON, ANTIPROTON);

      const FinalState fs(-2.0, 2.0);
      addProjection(fs, "FS");
      addProjection(FastJets(fs), "Jets"); 
      addProjection(TotalVisibleMomentum(fs), "CalMET");
      addProjection(PVertex(), "PV");
      // Veto (anti)neutrinos, and muons with pT above 1.0 GeV
      VetoedFinalState vfs(fs);
      vfs
        .addVetoPairId(NU_E)
        .addVetoPairId(NU_MU)
        .addVetoPairId(NU_TAU)
        .addVetoDetail(MUON, 1.0*GeV, MAXDOUBLE);
      addProjection(vfs, "VFS");
      addProjection(JetShape(vfs, _jetaxes, 0.0, 0.7, 0.1, 0.3, ENERGY), "JetShape");

      // Specify pT bins and initialise weight entries
      /// @todo Get these numbers from bundled data files
      _pTbins.resize(19);
      double ptbins[] = { 37.0, 45.0, 55.0, 63.0, 73.0, 84.0, 97.0, 112.0, 128.0, 148.0, 
                          166.0, 186.0, 208.0, 229.0, 250.0, 277.0, 304.0, 340.0, 380.0 };
      for (size_t i = 0; i <= 18; ++i) {
        _pTbins[i]  =  ptbins[i];
        _ShapeWeights[i%18] = 0.0;
      }
    }
      

  public:

    /// Factory method
    static Analysis* create() { 
      return new CDF_2005_S6217184(); 
    }


    /// @name Publication metadata
    //@{
    /// Get SPIRES ID code.
    string getSpiresId() const {
      return "6217184";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "CDF Run II jet shape paper";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2005";
    }
    /// Journal, and preprint references.
    virtual vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("hep-ex/0505013");
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

    /// @name Analysis cuts
    //@{
    /// Cut on primary vertex z-position (\f$ z(\text{PV}) < 60 \text{cm} \f$)
    const double _pvzmax;
    //@}

    /// @name Analysis data
    //@{
    vector<FourMomentum> _jetaxes;

    double _Rjet;

    double _ShapeWeights[18];

    /// \f$p_\perp\f$ bins to be distinguished during analysis
    vector<double> _pTbins;
    //@}


    //@{
    /// Histograms
    AIDA::IProfile1D* _profhistRho_pT[18];
    AIDA::IProfile1D* _profhistPsi_pT[18];
    AIDA::IProfile1D* _profhistPsi;
    //@}


  private:
    /// Hide the assignment operator
    CDF_2005_S6217184& operator=(const CDF_2005_S6217184&);

  };

}

#endif
