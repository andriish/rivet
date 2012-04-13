// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include "LWH/AIManagedObject.h"
using namespace AIDA;

namespace Rivet {


  /// @brief SLD flavour-dependent fragmentation paper
  /// @author Peter Richardson
  class SLD_2004_S5693039 : public Analysis {
  public:

    /// Constructor
    SLD_2004_S5693039() : Analysis("SLD_2004_S5693039"),
	_weightedTotalChargedPartNumLight(0.),
	_weightedTotalChargedPartNumCharm(0.),
	_weightedTotalChargedPartNumBottom(0.),
	_weightLight(0.),_weightCharm(0.),_weightBottom(0.)
    {}

    /// @name Analysis methods
    //@{

    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 2 charged FS particles
      const FinalState& fs = applyProjection<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed ncharged cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");
      // Get event weight for histo filling
      const double weight = e.weight();
   
      // Get beams and average beam momentum
      const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.momentum().vector3().mod() +
                                   beams.second.momentum().vector3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      int flavour = 0;
      const InitialQuarks& iqf = applyProjection<InitialQuarks>(e, "IQF");
   
      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = abs( iqf.particles().front().pdgId() );
      }
      else {
        map<int, double> quarkmap;
        foreach (const Particle& p, iqf.particles()) {
          if (quarkmap[p.pdgId()] < p.momentum().E()) {
            quarkmap[p.pdgId()] = p.momentum().E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }
      // total multiplicities
      switch (flavour) {
      case 1: case 2: case 3:
	_weightLight  += weight;
        _weightedTotalChargedPartNumLight  += numParticles * weight;
        break;
      case 4:
	_weightCharm  += weight;
        _weightedTotalChargedPartNumCharm  += numParticles * weight;
        break;
      case 5:
	_weightBottom += weight;
        _weightedTotalChargedPartNumBottom += numParticles * weight;
        break;
      }
      // spectra and individuall multiplicities
      foreach (const Particle& p, fs.particles()) {
	double pcm = p.momentum().vector3().mod();
	const double xp = pcm/meanBeamMom;

	_histPCharged ->fill(pcm     , weight);
	// all charged
	switch (flavour) {
	case DQUARK: case UQUARK: case SQUARK:
	  _histXpChargedL->fill(xp, weight);
	  break;
	case CQUARK:
	  _histXpChargedC->fill(xp, weight);
	  break;
	case BQUARK:
	  _histXpChargedB->fill(xp, weight);
	  break;
	}

	int id = abs(p.pdgId());
	// charged pions
	if(id==PIPLUS) {
	  _histXpPiPlus->fill(xp, weight);
	  _histXpPiPlusTotal->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _histXpPiPlusL->fill(xp, weight);
	    _multPiPlusL->fill(sqrtS(), weight);
	    break;
	  case CQUARK:
	    _histXpPiPlusC->fill(xp, weight);
	    _multPiPlusC->fill(sqrtS(), weight);
	    break;
	  case BQUARK:
	    _histXpPiPlusB->fill(xp, weight);
	    _multPiPlusB->fill(sqrtS(), weight);
	    break;
	  }
	}
	else if(id==KPLUS) {
	  _histXpKPlus->fill(xp, weight);
	  _histXpKPlusTotal->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _histXpKPlusL->fill(xp, weight);
	    _multKPlusL->fill(sqrtS(), weight);
	    break;
	  case CQUARK:
	    _histXpKPlusC->fill(xp, weight);
	    _multKPlusC->fill(sqrtS(), weight);
	    break;
	  case BQUARK:
	    _histXpKPlusB->fill(xp, weight);
	    _multKPlusB->fill(sqrtS(), weight);
	    break;
	  }
	}
	else if(id==PROTON) {
	  _histXpProton->fill(xp, weight);
	  _histXpProtonTotal->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _histXpProtonL->fill(xp, weight);
	    _multProtonL->fill(sqrtS(), weight);
	    break;
	  case CQUARK:
	    _histXpProtonC->fill(xp, weight);
	    _multProtonC->fill(sqrtS(), weight);
	    break;
	  case BQUARK:
	    _histXpProtonB->fill(xp, weight);
	    _multProtonB->fill(sqrtS(), weight);
	    break;
	  }
	}
      }
    }

    void init() {
      // Projections
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "FS");
      addProjection(InitialQuarks(), "IQF"); 
      _histPCharged   = bookHistogram1D( 1, 1, 1);
      _histXpPiPlus   = bookHistogram1D( 2, 1, 2);
      _histXpKPlus    = bookHistogram1D( 3, 1, 2);
      _histXpProton   = bookHistogram1D( 4, 1, 2);
      _histXpPiPlusTotal = bookHistogram1D( 2, 2, 2);
      _histXpKPlusTotal  = bookHistogram1D( 3, 2, 2);
      _histXpProtonTotal = bookHistogram1D( 4, 2, 2);
      _histXpPiPlusL  = bookHistogram1D( 5, 1, 1);
      _histXpPiPlusC  = bookHistogram1D( 5, 1, 2);
      _histXpPiPlusB  = bookHistogram1D( 5, 1, 3);
      _histXpKPlusL   = bookHistogram1D( 6, 1, 1);
      _histXpKPlusC   = bookHistogram1D( 6, 1, 2);
      _histXpKPlusB   = bookHistogram1D( 6, 1, 3);
      _histXpProtonL  = bookHistogram1D( 7, 1, 1);
      _histXpProtonC  = bookHistogram1D( 7, 1, 2);
      _histXpProtonB  = bookHistogram1D( 7, 1, 3);
      _histXpChargedL = bookHistogram1D( 8, 1, 1);
      _histXpChargedC = bookHistogram1D( 8, 1, 2);
      _histXpChargedB = bookHistogram1D( 8, 1, 3);

      _multPiPlusL  = bookHistogram1D( 5, 2, 1);
      _multPiPlusC  = bookHistogram1D( 5, 2, 2);
      _multPiPlusB  = bookHistogram1D( 5, 2, 3);
      _multKPlusL   = bookHistogram1D( 6, 2, 1);
      _multKPlusC   = bookHistogram1D( 6, 2, 2);
      _multKPlusB   = bookHistogram1D( 6, 2, 3);
      _multProtonL  = bookHistogram1D( 7, 2, 1);
      _multProtonC  = bookHistogram1D( 7, 2, 2);
      _multProtonB  = bookHistogram1D( 7, 2, 3);

    }


    /// Finalize
    void finalize() {
      // multiplicities
      // bottom
      const double avgNumPartsBottom = _weightedTotalChargedPartNumBottom / _weightBottom;
      AIDA::IDataPointSet * multB = bookDataPointSet(8, 2, 3);
      multB->point(0)->coordinate(1)->setValue(avgNumPartsBottom);
      // charm
      const double avgNumPartsCharm = _weightedTotalChargedPartNumCharm / _weightCharm;
      AIDA::IDataPointSet * multC = bookDataPointSet(8, 2, 2);
      multC->point(0)->coordinate(1)->setValue(avgNumPartsCharm);
      // light
      const double avgNumPartsLight = _weightedTotalChargedPartNumLight / _weightLight;
      AIDA::IDataPointSet * multL = bookDataPointSet(8, 2, 1);
      multL->point(0)->coordinate(1)->setValue(avgNumPartsLight);
      // charm-light
      AIDA::IDataPointSet * multD1 = bookDataPointSet(8, 3, 2);
      multD1->point(0)->coordinate(1)->setValue(avgNumPartsCharm -avgNumPartsLight);
      // bottom-light
      AIDA::IDataPointSet * multD2 = bookDataPointSet(8, 3, 3);
      multD2->point(0)->coordinate(1)->setValue(avgNumPartsBottom-avgNumPartsLight);
      // histograms 
      scale(_histPCharged   ,1./sumOfWeights());
      scale(_histXpPiPlus   ,1./sumOfWeights());
      scale(_histXpKPlus    ,1./sumOfWeights());
      scale(_histXpProton   ,1./sumOfWeights());
      scale(_histXpPiPlusTotal ,1./sumOfWeights());
      scale(_histXpKPlusTotal  ,1./sumOfWeights());
      scale(_histXpProtonTotal ,1./sumOfWeights());
      scale(_histXpPiPlusL  ,1./_weightLight);
      scale(_histXpPiPlusC  ,1./_weightCharm);
      scale(_histXpPiPlusB  ,1./_weightBottom);
      scale(_histXpKPlusL   ,1./_weightLight);
      scale(_histXpKPlusC   ,1./_weightCharm);
      scale(_histXpKPlusB   ,1./_weightBottom);
      scale(_histXpProtonL  ,1./_weightLight);
      scale(_histXpProtonC  ,1./_weightCharm);
      scale(_histXpProtonB  ,1./_weightBottom);

      scale(_histXpChargedL ,1./_weightLight);
      scale(_histXpChargedC ,1./_weightCharm);
      scale(_histXpChargedB ,1./_weightBottom);

      scale(_multPiPlusL   ,1./_weightLight);
      scale(_multPiPlusC   ,1./_weightCharm);
      scale(_multPiPlusB   ,1./_weightBottom);
      scale(_multKPlusL    ,1./_weightLight);
      scale(_multKPlusC    ,1./_weightCharm);
      scale(_multKPlusB    ,1./_weightBottom);
      scale(_multProtonL   ,1./_weightLight);
      scale(_multProtonC   ,1./_weightCharm);
      scale(_multProtonB   ,1./_weightBottom);

    }
    //@}

  private:


    /// @name Multiplicities
    //@{
    double _weightedTotalChargedPartNumLight;
    double _weightedTotalChargedPartNumCharm;
    double _weightedTotalChargedPartNumBottom;
    //@}

    /// @name Weights
    //@{
    double _weightLight;
    double _weightCharm;
    double _weightBottom;
    //@}

    // Histograms
    //@{
    AIDA::IHistogram1D *_histPCharged  ;
    AIDA::IHistogram1D *_histXpPiPlus  ;
    AIDA::IHistogram1D *_histXpKPlus   ;
    AIDA::IHistogram1D *_histXpProton  ;
    AIDA::IHistogram1D *_histXpPiPlusTotal;
    AIDA::IHistogram1D *_histXpKPlusTotal ;
    AIDA::IHistogram1D *_histXpProtonTotal;
    AIDA::IHistogram1D *_histXpPiPlusL ;
    AIDA::IHistogram1D *_histXpPiPlusC ;
    AIDA::IHistogram1D *_histXpPiPlusB ;
    AIDA::IHistogram1D *_histXpKPlusL  ;
    AIDA::IHistogram1D *_histXpKPlusC  ;
    AIDA::IHistogram1D *_histXpKPlusB  ;
    AIDA::IHistogram1D *_histXpProtonL ;
    AIDA::IHistogram1D *_histXpProtonC ;
    AIDA::IHistogram1D *_histXpProtonB ;
    AIDA::IHistogram1D *_histXpChargedL;
    AIDA::IHistogram1D *_histXpChargedC;
    AIDA::IHistogram1D *_histXpChargedB;
    AIDA::IHistogram1D *_multPiPlusL ;
    AIDA::IHistogram1D *_multPiPlusC ;
    AIDA::IHistogram1D *_multPiPlusB ;
    AIDA::IHistogram1D *_multKPlusL  ;
    AIDA::IHistogram1D *_multKPlusC  ;
    AIDA::IHistogram1D *_multKPlusB  ;
    AIDA::IHistogram1D *_multProtonL ;
    AIDA::IHistogram1D *_multProtonC ;
    AIDA::IHistogram1D *_multProtonB ;
    //@}
    
  };  
      
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SLD_2004_S5693039);

}
