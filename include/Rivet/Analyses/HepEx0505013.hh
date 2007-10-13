// -*- C++ -*-
#ifndef RIVET_HepEx0505013_HH
#define RIVET_HepEx0505013_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/KtJets.hh"
#ifdef HAVE_FASTJET	
#include "Rivet/Projections/FastJets.hh"
#endif
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/JetShape.hh"
#include "Rivet/RivetAIDA.fhh"
#include "Rivet/RivetCLHEP.hh"


namespace Rivet {

  /// Implementation of CDF RunII jet shape paper hep-ex/0505013 
  class HepEx0505013 : public Analysis {

  public:

    /// Default constructor
    inline HepEx0505013()
      // NB. eta in [-2.,2.] cut specified via FinalState constructor
      // NB. jetshape rmin=0.0, rmax=0.7, interval=0.1, r1minPsi=0.3
      : _fsproj(-2., 2.), _vfsproj(_fsproj), _jetsproj(_fsproj), 
        _calmetproj(_fsproj), _vertexproj(), 
        _jetshapeproj(_vfsproj, _jetaxes, 0.0, 0.7, 0.1, 0.3, ENERGY) 
    { 

      setBeams(PROTON, ANTIPROTON);


      // Add particle/antiparticle vetoing: 12=nu_e, 14=nu_mu, 16=nu_tau
      _vfsproj
        .addVetoPairId(12)
        .addVetoPairId(14)
        .addVetoPairId(16);
      
      // Veto muons (PDG code = 13) with pT above 1.0 GeV
      _vfsproj.addVetoDetail(13, 1.0, numeric_limits<double>::max());

      addProjection(_fsproj);
      addProjection(_vfsproj);
      addProjection(_jetsproj);
      addProjection(_calmetproj);
      addProjection(_vertexproj);
      addProjection(_jetshapeproj);


      _Rjet = 0.7;

      _pTbins.resize(19);
      _pTbins[0]  =  37.;
      _pTbins[1]  =  45.;
      _pTbins[2]  =  55.;
      _pTbins[3]  =  63.;
      _pTbins[4]  =  73.;
      _pTbins[5]  =  84.;
      _pTbins[6]  =  97.;
      _pTbins[7]  = 112.;
      _pTbins[8]  = 128.;
      _pTbins[9]  = 148.;
      _pTbins[10] = 166.;
      _pTbins[11] = 186.;
      _pTbins[12] = 208.;
      _pTbins[13] = 229.;
      _pTbins[14] = 250.;
      _pTbins[15] = 277.;
      _pTbins[16] = 304.;
      _pTbins[17] = 340.;
      _pTbins[18] = 380.;

      for (int i=0; i<18; ++i)
	_ShapeWeights[i] = 0.;

    }
      

  public:

    /// Factory method
    static Analysis* create() { return new HepEx0505013(); }

    /// Return the name of this analysis.
    inline string getName() const {
      return "HepEx0505013";
    }

  public:

    void init();
    
    void analyze(const Event & event);
    
    void finalize();

  private:

    /// The final state projector used by this analysis.
    FinalState _fsproj;
    
    ///The vetoed final state projector needed by the jet algorithm
    VetoedFinalState _vfsproj; 
    
    /// The D0ILConeJets projector used by this analysis.
    ////D0ILConeJets _jetsproj;
    /// The KtJets projector alternatively usable by this analysis.
    ////KtJets _jetsproj;
    #ifdef HAVE_FASTJET
    /// The FastJets projector alternatively usable by this analysis.
    FastJets _jetsproj;
    #else
    /// The D0ILConeJets projector used by this analysis.
    D0ILConeJets _jetsproj;
    /// The KtJets projector alternatively usable by this analysis.
    ////KtJets _jetsproj;
    #endif
    
    /// The Calorimeter Missing Et projector
    TotalVisibleMomentum _calmetproj;

    /// The Primary Vertex projector
    PVertex _vertexproj;

    /// The Calorimeter Missing Et projector
    JetShape _jetshapeproj;
    

  private:

    vector<LorentzVector> _jetaxes;

    double _Rjet;

    double _ShapeWeights[18];

    /// pT bins to be distiguished during analysis
    vector<double> _pTbins;


    /// Hide the assignment operator
    HepEx0505013& operator=(const HepEx0505013&);

    //@{
    /// Histograms
    AIDA::IProfile1D* _profhistRho_pT[18];
    AIDA::IProfile1D* _profhistPsi_pT[18];
    AIDA::IProfile1D* _profhistPsi;

    //@}

  };

}

#endif
