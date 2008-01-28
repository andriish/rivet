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
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  
  /// This CDF analysis provides \f$ p_T \f$ and \f$ \eta \f$ distributions 
  /// of jets in Z + (b) jet production, before and after tagging.
  ///
  /// Analysis cut values:
  ///  - _pvzmax: cut on primary vertex \f$ z \f$ position (\f$ z(\text{PV}) < 60 \text{cm} \f$);
  ///  - _met: cut on missing transverse energy (\f$ MET < 25 \text{GeV} \f$)
  ///  - _set: cut on scalar sum of transverse energy (\f$ SET < 150 \text{GeV} \f$)
  ///  - _lepPtmin: cut on leptons minimum transverse momentum (\f$ pT(\text{lep.}) > 10 \text{GeV} \f$)
  ///  - _eletamax: cut on max. pseudorapidity of electrons/calorimeter (\f$ |eta(el.)| < 3.6 \f$) 
  ///  - _muetamax: cut on max. pseudorapidity of muons/muon chambers (\f$ |eta(el.)| < 1.5 \f$) 
  ///  - _mllmin: cut on min invariant dilepton mass (\f$ m(ll) > 66 \text{GeV} \f$)
  ///  - _mllmax: cut on max invariant dilepton mass (\f$ m(ll) <116 \text{GeV} \f$)
  ///  - _triglepPtmin: cut on trigger lepton min. transverse momentum (\f$ pT > 18 \text{GeV} \f$)
  ///  - _triglepetamax: cut on trigger lepton max pseudorapidity (\f$ |eta| < 1.1 \f$)
  ///  - _Rlepisol: cut on lepton isolation radius (\f$ R < 0.4 \f$)
  ///  - hadronic calorimeter fraction of charged hadrons, needed for cuts on hadr. + el. mag. fractions
  ///  - hadronic/el.magn. calo. fraction cut const (fhfemconst = 0.055)
  ///  - hadronic/el.magn. calo. fraction cut slope (fhfemslope = 0.00045)
  ///  - muon Energy cut case separation (muEsep = 100 \text{GeV})
  ///  - muon min. E hadr.isol. threshold (muEhMin = 6 \text{GeV})
  ///  - muon min. E el.mag. isol. threshold (muEemMin = 2 \text{GeV})
  ///  - muon had. Energy isol. cut slope (muEhslope = 0.0280)
  ///  - muon el.mag. Energy isol. cut slope (muEemslope = 0.0115 \f$)
  ///  - jets transverse momentum cut (\f$ p_\text{T,jet} > 20 \text{GeV} \f$)
  ///  - jets pseudorapidity cut (\f$ \eta_\text{jet} < 1.5 \f$)
  ///  - _vtxjetRmax: Max. distance in \f$ (\eta, \phi) \f$ space between vertex vis. momentum and jet to be probed (\f$ R_\text{vtx-jet} < 0.7 \f$) 
  ///  - _trketamax: Tracker geometrical acceptance (\f$ \eta < 2.0 \f$)
  ///  - _ipres: Impact Parameter resolution (\f$ \Delta{\text{IP}} = 34e-3 \text{mm} \f$), including beam spot
  ///  - _dlsmin: cut on Decay Length Significance (\f$ l/\Delta{l} > 7.5 \f$)
  ///  - _dlsres: Decay Length Significance resolution (assumed to be (\f$ \Delta{l} = 34e-3 \text{mm} \f$))

  class CDF_2006_S6653332 : public Analysis {

  public:

    /// Default constructor
    CDF_2006_S6653332()
      : _fsproj(-3.6, 3.6), _vfsproj(_fsproj), _jetsproj(_vfsproj),
        _chfsproj(_vfsproj),
        _calmetproj(_vfsproj), _chleproj(_vfsproj), _pvtxproj(),
        _pvzmax(600*mm), _metmax(25*GeV), _setmax(150*GeV), _lepPtmin(10*GeV),
        _eletamax(3.5), _muetamax(1.5), _mllmin(66*GeV), _mllmax(116*GeV),
        _triglepPtmin(18*GeV), _trigeletamax(1.1), _trigmuetamax(1.0),
        _Rlepisol(0.4), _fh(0.8), _fhfemconst(0.055), _fhfemslope(0.00045),
        _muEsep(100*GeV), _muEhMin(6*GeV), _muEemMin(2*GeV), _muEhslope(0.0280),
        _muEemslope(0.0115), _jetPtmin(20*GeV), _jetetamax(1.5),
        _vtxjetRmax(0.7), _trketamax(2.0), _ipres(34e-3*mm), _dlsmin(7.5), _dlsres(34e-3*mm),
        _svtxproj(_pvtxproj, _chfsproj, _jetaxes, _vtxjetRmax, _trketamax, _ipres, _dlsmin, _dlsres)
    {

      setBeams(PROTON, ANTIPROTON);

      // Add particle/antiparticle vetoing: 12=nu_e, 14=nu_mu, 16=nu_tau
      _vfsproj
        .addVetoPairId(12)
        .addVetoPairId(14)
        .addVetoPairId(16);

      // Veto muons (PDG code = 13) with \f$ p_T \f$ above 1.0 GeV
      _vfsproj.addVetoDetail(13, 1.0, numeric_limits<double>::max());

      addProjection(_fsproj);
      addProjection(_vfsproj);
      addProjection(_chfsproj);
      addProjection(_jetsproj);
      addProjection(_calmetproj);
      addProjection(_chleproj);
      addProjection(_pvtxproj);
      addProjection(_svtxproj);
    }


    /// Factory method.
    static Analysis* create() { 
      return new CDF_2006_S6653332(); 
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "6653332";
    }

    /// Get a description of the analysis.
    string getDescription() const {
      return "";
    }

    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "CDF";
    }

    /// When published (preprint year according to SPIRES).
    string getYear() const {
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

    /// The final state projector.
    FinalState _fsproj;

    /// The visible final state projector.
    VetoedFinalState _vfsproj;

    /// The FastJets projector used by this analysis.
    FastJets _jetsproj;

    /// The charged final state projector.
    ChargedFinalState _chfsproj;

    /// The Calorimeter MET projector.
    TotalVisibleMomentum _calmetproj;

    /// The charged leptons projector.
    ChargedLeptons _chleproj;

    /// Projection to find the primary vertex.
    PVertex _pvtxproj;

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


  private:

    /// Set of vectors defining the jet axes (for detached vertex finding).
    vector<FourMomentum> _jetaxes;

    /// Projection to find secondary vertices.
    /// Needs to be initialized after constants on which constructor depends
    SVertex _svtxproj;
    
    /// Hide the assignment operator
    CDF_2006_S6653332& operator=(const CDF_2006_S6653332&);

    //@{
    /// Histograms
    AIDA::IHistogram1D* _histJetsPt;
    AIDA::IHistogram1D* _histJetsEta;
    AIDA::IHistogram1D* _histbJetsPt;
    AIDA::IHistogram1D* _histbJetsEta;
    //@}

  };

}

#endif
