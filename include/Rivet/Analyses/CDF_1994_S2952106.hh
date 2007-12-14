// -*- C++ -*-
#ifndef RIVET_CDF_1994_S2952106_HH
#define RIVET_CDF_1994_S2952106_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/RivetAIDA.fhh"



namespace Rivet {  

  /// Analysis based on the CDF Run I color coherence analysis described in PRD50, 9, 5562 (1994).
  /// >= three jet events are selected,
  /// Et distributions of the leading three pT jets are obtained.
  /// DeltaR between 2nd and 3rd leading jets in pT and pseudorapidiy of the 3rd leading jet
  /// are plotted. alpha = dH/dPhi is plotted, where dH is the pseudorapidity difference
  /// between the 2nd and 3rd leading jet and dPhi the azimuthal angle difference of these.
  /// Since the data has not been corrected to particle final state, a bin by bin correction is 
  /// applied, based on the distributions with ideal and CDF simulation as given in the publication.
  /// Analysis cut values:
  /// _pvzmax: cunt on primary vertex z position (z(PV) < 60 cm)
  /// _leadJetPt, _3rdJetPt: Min. Pt of the leading and 3rd leading jets
  /// _etamax: Max. pseudorapidity range of 2nd and 3rd leading jets
  /// _phimin: Delta phi (azimuthal angle) requirement (transverse back to back'ness)
  /// _metsetmax: MET over sqrt(Scalar ET) cut requirement




  /// @author Lars Sonnenschein
  class CDF_1994_S2952106 : public Analysis {

  public:

    /// Default constructor.
    inline CDF_1994_S2952106()
      // NB. eta in [-4.2, 4.2] cut specified via FinalState constructor, CDF CAL acceptance
      : _fsproj(-4.2, 4.2), _vfsproj(_fsproj), 

	//_conejetsproj(fastjet::kt_algorithm, fastjet::Et_scheme, 0.7, _fsproj),

	//_jetalgotype(FastJets::SISCone), _Rpar(0.7),
	//_conejetsproj(_jetalgotype, _Rpar, _fsproj), 

	_jetalgotype(FastJets::CDFJetClu), _Rpar(0.7),
	_conejetsproj(_jetalgotype, _Rpar, _fsproj), 

	_calmetproj(_fsproj), _vertexproj(),

	/// @todo z- value assumed to be in mm, PYTHIA convention: dangerous!
	_pvzmax(600.), _leadJetPt(100.), _etamax(0.7), _phimin(1./9.*PI/2.), 
	_metsetmax(6.), _3rdJetPt(10.)

    { 
      
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);


      // Add particle/antiparticle vetoing: 12=nu_e, 14=nu_mu, 16=nu_tau
      _vfsproj
        .addVetoPairId(12)
        .addVetoPairId(14)
        .addVetoPairId(16);
      
      // Veto muons (PDG code = 13) with pT above 1.0 GeV
      _vfsproj.addVetoDetail(13, 1.0, numeric_limits<double>::max());


      addProjection(_fsproj);
      addProjection(_vfsproj);
      addProjection(_conejetsproj);
      addProjection(_calmetproj);
      addProjection(_vertexproj);
      
    }

    /// Factory method
    static Analysis* create() { return new CDF_1994_S2952106(); }

    /// Return the name of this analysis.
    inline string getName() const {
      return "CDF_1994_S2952106";
    }

  public:

    void init();
    
    void analyze(const Event& event);
    
    void finalize();

  private:

    /// The final state projector used by this analysis.
    FinalState _fsproj;

    ///The vetoed final state projector needed by the jet algorithm
    VetoedFinalState _vfsproj; 


    /// Parameters used in the CDF Run I cone jet algorithm
    //fastjet::JetFinder _jetalgo; //=plugin_algorithm
    //fastjet::RecombinationScheme _recomb; //external_scheme
    FastJets::ExtJetType _jetalgotype;

    const double _Rpar; // =0.7
    ///The Fastjet plugin CDF Run I cone jet algorithm
    FastJets _conejetsproj;

    /// The D0ILConeJets projector used by this analysis.
    //D0ILConeJets _conejetsproj;

    /// The Calorimeter Missing Et projector
    TotalVisibleMomentum _calmetproj;

    /// The Primary Vertex projector
    PVertex _vertexproj;

    /// Counter for the number of events analysed
    double _eventsTried;

    // Analysis cuts
    ///Cut on primary vertex z-position (z(PV)<60 cm)
    const double _pvzmax;
    /// Min. Pt of the leading and 3rd leading jets
    const double _leadJetPt;
    const double _3rdJetPt;
    /// Max. pseudorapidity range of 2nd and 3rd leading jets
    const double _etamax;
    /// Delta phi (azimuthal angle) requirement (transverse back to back'ness)
    const double _phimin;
    /// MET over sqrt(Scalar ET) cut requirement
    const double _metsetmax;

    /// Hide the assignment operator
    CDF_1994_S2952106& operator=(const CDF_1994_S2952106& x);

    /// @name Histograms
    //@{
    // AIDA::IHistogram2D* _histHvsDphi;
    // AIDA::IHistogram2D* _histRvsAlpha;
    AIDA::IHistogram1D* _histJet1Et;
    AIDA::IHistogram1D* _histJet2Et;
    AIDA::IHistogram1D* _histR23;
    AIDA::IHistogram1D* _histJet3eta;
    AIDA::IHistogram1D* _histAlpha;
    // AIDA::IHistogram1D* _histAlphaMCvsDat;
    AIDA::IHistogram1D* _histAlpaIdeal;
    AIDA::IHistogram1D* _histAlpaCDF;
    AIDA::IHistogram1D* _histR23Ideal;
    AIDA::IHistogram1D* _histR23CDF;
    AIDA::IHistogram1D* _histJet3etaIdeal;
    AIDA::IHistogram1D* _histJet3etaCDF;
    //@}

  };

}

#endif
