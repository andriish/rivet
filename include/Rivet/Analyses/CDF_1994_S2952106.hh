// -*- C++ -*-
#ifndef RIVET_CDF_1994_S2952106_HH
#define RIVET_CDF_1994_S2952106_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {  

  /// Analysis based on the CDF Run I color coherence analysis described in
  /// PRD50,9,5562 (1994). Events with >= 3 jets are selected and \f$ E_T \f$
  /// distributions of the leading three \f$ p_T \f$ jets are obtained.  \f$
  /// \Delta{R} \f$ between 2nd and 3rd leading jets in \f$ p_T \f$ and
  /// pseudorapidity of the 3rd leading jet are plotted. \f$ \alpha =
  /// \d{\eta}/\d{\phi} \f$ is plotted, where \f$ \d{\eta} \f$ is the
  /// pseudorapidity difference between the 2nd and 3rd leading jet and dPhi the
  /// azimuthal angle difference of these.  Since the data has not been
  /// corrected to particle final state, a bin by bin correction is applied,
  /// based on the distributions with ideal and CDF simulation as given in the
  /// publication.
  ///
  /// Analysis cut values: 
  ///  - _pvzmax: cut on primary vertex \f$ z \f$ position (\f$ z(\text{PV}) < 60 \text{cm} \f$);
  ///  - _leadJetPt, _3rdJetPt: Min. \f$ p_T \f$ of the leading and 3rd leading jets;
  ///  - _etamax: Max. pseudorapidity range of 2nd and 3rd leading jets;
  ///  - _phimin: \f$ \Delta{\phi} \f$ (azimuthal angle) requirement (transverse back to back'ness);
  ///  - _metsetmax: MET over \f$ \sqrt{\text{Scalar }E_T} \f$ cut requirement.
  ///
  /// @author Lars Sonnenschein
  ///
  class CDF_1994_S2952106 : public Analysis {

  public:

    /// Constructor.
    /// NB. eta in [-4.2, 4.2] cut specified via FinalState constructor, CDF CAL acceptance.
	/// @todo Set units on _metsetmax
    CDF_1994_S2952106()
      : _fsproj(-4.2, 4.2), _vfsproj(_fsproj),
        _conejetsproj(_fsproj, FastJets::CDFJETCLU, 0.7),
        _calmetproj(_fsproj), _vertexproj(),
        _pvzmax(600*mm), _leadJetPt(100*GeV), _3rdJetPt(10*GeV),
        _etamax(0.7), _phimin(PI/18.0), _metsetmax(6.0)
    {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);

      // Add particle/antiparticle vetoing: 12=nu_e, 14=nu_mu, 16=nu_tau
      /// @todo Use ParticleName enum
      _vfsproj
        .addVetoPairId(12)
        .addVetoPairId(14)
        .addVetoPairId(16);
      
      // Veto muons with pT above 1.0 GeV
      _vfsproj.addVetoDetail(MUON, 1.0*GeV, MAXDOUBLE);

      addProjection(_fsproj);
      addProjection(_vfsproj);
      addProjection(_conejetsproj);
      addProjection(_calmetproj);
      addProjection(_vertexproj);
    }


    /// Factory method
    static Analysis* create() { return new CDF_1994_S2952106(); }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "2952106";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "CDF Run I color coherence analysis.";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "1994";
    }
    /// Publication references.
    vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("Phys.Rev.D50,9,5562 (1994)");
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

    /// The final state projector used by this analysis.
    FinalState _fsproj;

    /// The vetoed final state projector needed by the jet algorithm.
    VetoedFinalState _vfsproj; 

    ///The cone jet algorithm.
    FastJets _conejetsproj;
    //D0ILConeJets _conejetsproj;

    /// The Calorimeter Missing Et projector
    TotalVisibleMomentum _calmetproj;

    /// The Primary Vertex projector
    PVertex _vertexproj;

    /// Counter for the number of events analysed
    double _eventsTried;

    /// @name Analysis cuts
    //@{
    ///Cut on primary vertex z-position (z(PV) < 60 cm)
    const double _pvzmax;

    /// Min \f$ p_T \f$ of the leading and 3rd leading jets.
    //@{
    const double _leadJetPt;
    const double _3rdJetPt;
    //@}

    /// Max pseudorapidity range of 2nd and 3rd leading jets.
    const double _etamax;

    /// Delta phi (azimuthal angle) requirement (transverse back to back'ness).
    const double _phimin;

    /// MET over sqrt(scalar \f$ E_T \f$) cut requirement.
    const double _metsetmax;
    //@}

  private:
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
