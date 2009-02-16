// -*- C++ -*-
#ifndef RIVET_CDF_2005_S6217184_HH
#define RIVET_CDF_2005_S6217184_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /* @brief CDF Run II jet shape analysis
   * @author Lars Sonnenschein
   *
   *
   * @par Run conditions
   *
   * @arg \f$ \sqrt{s} = \f$ 1960 GeV
   * @arg Run with generic QCD events.
   * @arg \f$ \eta \in [-2,2] \f$ cut used on final state.
   * @arg Jet shape \f$ r_\text{min} = 0.0 \f$, \f$ r_\text{max} = 0.7 \f$,
   * @arg radial interval = 0.1, r1minPsi = 0.3.
   *
   */	
  class CDF_2005_S6217184 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    // Constructor
    /// Cuts on final state are \f$ \eta \in [-2,2] \f$.
    /// Jet shape \f$ r_\text{min} = 0.0 \f$, \f$ r_\text{max} = 0.7 \f$, 
    /// interval = 0.1, r1minPsi = 0.3.
    CDF_2005_S6217184();

    /// Factory method
    static Analysis* create() { 
      return new CDF_2005_S6217184(); 
    }

    //@}

    /// @name Publication metadata
    //@{
    /// Get SPIRES ID code.
    string spiresId() const {
      return "6217184";
    }
    /// A short description of the analysis.
    string summary() const {
      return "CDF Run II jet shape analysis";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2005";
    }
    /// Journal, and preprint references.
    virtual vector<string> references() const {
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


    /// @name Histogram collections
    //@{
    AIDA::IProfile1D* _profhistRho_pT[18];
    AIDA::IProfile1D* _profhistPsi_pT[18];
    AIDA::IProfile1D* _profhistPsi;
    //@}

  };


}

#endif
