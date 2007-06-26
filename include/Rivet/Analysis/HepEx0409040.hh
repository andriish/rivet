// -*- C++ -*-
#ifndef RIVET_HepEx0409040_H
#define RIVET_HepEx0409040_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {  

  /// Analysis based on the D0 Run II jet analysis described in hep-ex/0409040.
  /// @author Lars Sonnenschein
  class HepEx0409040 : public Analysis {

  public:

    /// Default constructor.
    inline HepEx0409040()
      // NB. eta in [-3,3] cut specified via FinalState constructor
      : _fsproj(-3.0, 3.0), _vfsproj(_fsproj), _calmetproj(_fsproj), _conejetsproj(_fsproj), _vertexproj()
    { 
      setBeams(PROTON, ANTIPROTON);
      addProjection(_fsproj);
      addProjection(_vertexproj);
      addProjection(_vfsproj);
      addProjection(_conejetsproj);

      // Add particle/antiparticle vetoing: 12=nu_e, 14=nu_mu, 16=nu_tau
      _vfsproj
        .addVetoPairId(12)
        .addVetoPairId(14)
        .addVetoPairId(16);
      
      // Veto muons (PDG code = 13) with pT above 1.0 GeV
      _vfsproj.addVetoDetail(13, 1.0, numeric_limits<double>::max());

      // Don't put neutrinos or low pT muons into the cal missing ET.
      /// @todo Spot the memory leak!
      _calmetproj = TotalVisibleMomentum(_vfsproj);
      addProjection(_calmetproj);
   }


    /// Return the name of this analysis.
    inline string getName() const {
      return "HepEx0409040";
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

    /// The Calorimeter Missing Et projector
    TotalVisibleMomentum _calmetproj;

    /// The D0ILConeJets projector used by this analysis.
    D0ILConeJets _conejetsproj;

    /// The Primary Vertex projector
    PVertex _vertexproj;

    /// Hide the assignment operator
    HepEx0409040& operator=(const HepEx0409040& x);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histJetAzimuth_pTmax75_100;
    AIDA::IHistogram1D* _histJetAzimuth_pTmax100_130;
    AIDA::IHistogram1D* _histJetAzimuth_pTmax130_180;
    AIDA::IHistogram1D* _histJetAzimuth_pTmax180_;
    //@}

  };

}

#endif
