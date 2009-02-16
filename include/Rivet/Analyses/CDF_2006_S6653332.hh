// -*- C++ -*-
#ifndef RIVET_CDF_2006_S6653332_HH
#define RIVET_CDF_2006_S6653332_HH

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


  /* @brief CDF Run II analysis: jet \f$ p_T \f$ and \f$ \eta \f$ distributions in Z + (b) jet production
   * @author Lars Sonnenschein
   *
   * This CDF analysis provides \f$ p_T \f$ and \f$ \eta \f$ distributions of
   * jets in Z + (b) jet production, before and after tagging.
   *
   * @par Run conditions
   *
   * @arg \f$ \sqrt{s} = \f$ 1960 GeV
   * @arg Run with Z+jets events.
   * @arg Jets transverse momentum cut ( \f$ p_\text{T,jet} > 20 \text{GeV} \f$ );
   * @arg Cut on leptons minimum transverse momentum ( \f$ pT(\text{lep.}) > 10 \text{GeV} \f$ );
   *
   */
  class CDF_2006_S6653332 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{
    
    /// Constructor
    CDF_2006_S6653332();

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

    /// A short description of the analysis.
    string description() const {
      return "pT and eta distributions of jets in Z + (b)jet production, before and after tagging.";
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }

    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2006";
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    /// @name Analysis cuts
    //@{
    /// Cut on primary vertex z-position (\f$ z(\text{PV}) < 60 \text{cm} \f$)
    const double _pvzmax;

    /// Cut on missing transverse energy (\f$ \text{MET} < 25 \text{GeV} \f$)
    const double _metmax;

    /// Cut on scalar sum of transverse energy (\f$ \text{SET} < 150 \text{GeV} \f$)
    const double _setmax;

    /// Cut on leptons minimum transverse momentum (\f$ p_T(\text{lep.}) > 10 \text{GeV} \f$) 
    const double _lepPtmin;

    /// Cut on max. pseudorapidity of electrons/calorimeter (\f$ |\eta(e)| < 3.6 \f$) 
    const double _eletamax;

    /// Cut on max. pseudorapidity of muons/muon chambers (\f$ |\eta(\mu)| < 1.5 \f$) 
    const double _muetamax;

    /// Cut on min invariant dilepton mass (\f$ m(ll) > 66 \text{GeV} \f$)
    const double _mllmin;

    /// Cut on max invariant dilepton mass (\f$ m(ll) < 116 \text{GeV} \f$)
    const double _mllmax;

    /// Cut on trigger lepton min transverse momentum (\f$ p_T > 18 \text{GeV} \f$)
    const double _triglepPtmin;

    /// Cut on trigger electron max pseudorapidity (\f$ |\eta(e)| < 1.1 \f$)
    const double _trigeletamax;

    /// Cut on trigger muon max pseudorapidity (\f$ |\eta(\mu)| < 1.0 \f$)
    const double _trigmuetamax;

    /// Cut on lepton isolation radius (\f$ R_\text{isol} < 0.4 \f$)
    const double _Rlepisol;

    /// Hadronic calorimeter fraction (fh = 0.8)
    const double _fh;

    /// Hadronic/el.magn. calo. fraction cut const (fhfemconst = 0.055)
    const double _fhfemconst;
    
    /// Hadronic/el.magn. calo. fraction cut slope (fhfemslope = 0.00045)
    const double _fhfemslope;

    /// Muon Energy cut case separation (muEsep = 100 GeV)
    const double _muEsep;

    /// Muon E hadr. isol. threshold (muEhMin = 6 GeV)
    const double _muEhMin;

    /// Muon E el.mag. isol. threshold (muEemMin = 2 GeV)
    const double _muEemMin;

    /// Muon had. Energy isol. cut slope (muEhslope = 0.0280)
    const double _muEhslope;

    /// Muon el.mag. Energy isol. cut slope (muEemslope = 0.0115)
    const double _muEemslope;
    
    /// Jets transverse momentum cut (jetPtmin = 20 GeV)
    const double _jetPtmin;

    /// Jets pseudorapidity cut (jetetamax = 1.5)
    const double _jetetamax;

    /// Max. distance (eta,phi) between vertex vis. momentum and jet to be probed (vtxjetRmax= 0.7) 
    const double _vtxjetRmax;

    /// Tracker geometrical acceptance (\f$ \eta < 2.0 \f$)
    const double _trketamax;

    /// Impact Parameter resolution (ipres = 34e-3mm), including beam spot
    const double _ipres;

    /// Cut on Decay Length Significance (dlsmin = 7.5)
    const double _dlsmin;

    /// Decay Length Significance resolution (assumed to be (dlsres = 34e-3mm))
    const double _dlsres;
    //@}


    /// Set of vectors defining the jet axes (for detached vertex finding).
    vector<FourMomentum> _jetaxes;

    
    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _histJetsPt;
    AIDA::IHistogram1D* _histJetsEta;
    AIDA::IHistogram1D* _histbJetsPt;
    AIDA::IHistogram1D* _histbJetsEta;
    //@}

  };


}

#endif
