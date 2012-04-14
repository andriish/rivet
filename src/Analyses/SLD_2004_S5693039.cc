// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include "Rivet/Projections/Thrust.hh"
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
      ParticleVector quarks;
      if (iqf.particles().size() == 2) {
        flavour = abs( iqf.particles().front().pdgId() );
	quarks = iqf.particles();
      }
      else {
        map<int, Particle > quarkmap;
        foreach (const Particle& p, iqf.particles()) {
          if (quarkmap.find(p.pdgId())==quarkmap.end())
	    quarkmap[p.pdgId()] = p;
	  else if (quarkmap[p.pdgId()].momentum().E() < p.momentum().E())
	    quarkmap[p.pdgId()] = p;
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
	  double energy(0.);
	  if(quarkmap.find( i)!=quarkmap.end())
	    energy += quarkmap[ i].momentum().E();
	  if(quarkmap.find(-i)!=quarkmap.end())
	    energy += quarkmap[-i].momentum().E();
          if (energy > maxenergy) flavour = i;
        }
	if(quarkmap.find( flavour)!=quarkmap.end())
	  quarks.push_back(quarkmap[ flavour]);
	if(quarkmap.find(-flavour)!=quarkmap.end())
	  quarks.push_back(quarkmap[-flavour]);
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
      // thrust axis for projections
      Vector3 axis = applyProjection<Thrust>(e, "Thrust").thrustAxis();
      double dot(0.);
      if(!quarks.empty()) {
	dot = quarks[0].momentum().vector3().dot(axis);
	if(quarks[0].pdgId()<0) dot *= -1.;
      }
      // spectra and individual multiplicities
      foreach (const Particle& p, fs.particles()) {
	double pcm = p.momentum().vector3().mod();
	const double xp = pcm/meanBeamMom;

	// if in quark or antiquark hemisphere
	bool quark = p.momentum().vector3().dot(axis)*dot>0.;

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
	    if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
	      _histRPiPlus->fill(xp, weight);
	    else
	      _histRPiMinus->fill(xp, weight);
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
	    if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
	      _histRKPlus->fill(xp, weight);
	    else
	      _histRKMinus->fill(xp, weight);
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
	    if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
	      _histRProton->fill(xp, weight);
	    else
	      _histRPBar  ->fill(xp, weight);
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
      addProjection(Thrust(FinalState()), "Thrust");
      // histograms 
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

      _histRPiPlus  = bookHistogram1D( 9, 1, 1);
      _histRPiMinus = bookHistogram1D( 9, 1, 2);
      _histRKPlus   = bookHistogram1D(10, 1, 1);
      _histRKMinus  = bookHistogram1D(10, 1, 2);
      _histRProton  = bookHistogram1D(11, 1, 1);
      _histRPBar    = bookHistogram1D(11, 1, 2);

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

      const string dir = histoDir();
      AIDA::IHistogram1D * num = histogramFactory().subtract(dir + "/n1",*_histRPiMinus,*_histRPiPlus);
      AIDA::IHistogram1D * den = histogramFactory().add     (dir + "/n2",*_histRPiMinus,*_histRPiPlus);
      AIDA::IDataPointSet* h   = histogramFactory().divide(dir +"/d09-x01-y03",*num,*den);
      scale(h,100.);
      histogramFactory().destroy(num);
      histogramFactory().destroy(den);
      num = histogramFactory().subtract(dir + "/n3",*_histRKMinus,*_histRKPlus);
      den = histogramFactory().add     (dir + "/n4",*_histRKMinus,*_histRKPlus);
      h   = histogramFactory().divide(dir +"/d10-x01-y03",*num,*den);
      scale(h,100.);
      histogramFactory().destroy(num);
      histogramFactory().destroy(den);
      num = histogramFactory().subtract(dir + "/n5",*_histRProton,*_histRPBar);
      den = histogramFactory().add     (dir + "/n6",*_histRProton,*_histRPBar);
      h   = histogramFactory().divide(dir +"/d11-x01-y03",*num,*den);
      scale(h,100.);
      histogramFactory().destroy(num);
      histogramFactory().destroy(den);

      // histograms 
      Analysis::scale(_histPCharged   ,1./sumOfWeights());
      Analysis::scale(_histXpPiPlus   ,1./sumOfWeights());
      Analysis::scale(_histXpKPlus    ,1./sumOfWeights());
      Analysis::scale(_histXpProton   ,1./sumOfWeights());
      Analysis::scale(_histXpPiPlusTotal ,1./sumOfWeights());
      Analysis::scale(_histXpKPlusTotal  ,1./sumOfWeights());
      Analysis::scale(_histXpProtonTotal ,1./sumOfWeights());
      Analysis::scale(_histXpPiPlusL  ,1./_weightLight);
      Analysis::scale(_histXpPiPlusC  ,1./_weightCharm);
      Analysis::scale(_histXpPiPlusB  ,1./_weightBottom);
      Analysis::scale(_histXpKPlusL   ,1./_weightLight);
      Analysis::scale(_histXpKPlusC   ,1./_weightCharm);
      Analysis::scale(_histXpKPlusB   ,1./_weightBottom);
      Analysis::scale(_histXpProtonL  ,1./_weightLight);
      Analysis::scale(_histXpProtonC  ,1./_weightCharm);
      Analysis::scale(_histXpProtonB  ,1./_weightBottom);

      Analysis::scale(_histXpChargedL ,1./_weightLight);
      Analysis::scale(_histXpChargedC ,1./_weightCharm);
      Analysis::scale(_histXpChargedB ,1./_weightBottom);

      Analysis::scale(_multPiPlusL   ,1./_weightLight);
      Analysis::scale(_multPiPlusC   ,1./_weightCharm);
      Analysis::scale(_multPiPlusB   ,1./_weightBottom);
      Analysis::scale(_multKPlusL    ,1./_weightLight);
      Analysis::scale(_multKPlusC    ,1./_weightCharm);
      Analysis::scale(_multKPlusB    ,1./_weightBottom);
      Analysis::scale(_multProtonL   ,1./_weightLight);
      Analysis::scale(_multProtonC   ,1./_weightCharm);
      Analysis::scale(_multProtonB   ,1./_weightBottom);

      // paper suggests this should be 0.5/weight but has to be 1.
      // to get normalisations right
      Analysis::scale(_histRPiPlus ,1./_weightLight);
      Analysis::scale(_histRPiMinus,1./_weightLight);
      Analysis::scale(_histRKPlus  ,1./_weightLight);
      Analysis::scale(_histRKMinus ,1./_weightLight);
      Analysis::scale(_histRProton ,1./_weightLight);
      Analysis::scale(_histRPBar   ,1./_weightLight);

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
    AIDA::IHistogram1D *_histRPiPlus ;
    AIDA::IHistogram1D *_histRPiMinus;
    AIDA::IHistogram1D *_histRKPlus  ;
    AIDA::IHistogram1D *_histRKMinus ;
    AIDA::IHistogram1D *_histRProton ;
    AIDA::IHistogram1D *_histRPBar   ;
    //@}
    
    void scale(AIDA::IDataPointSet*& histo, double scale) {
      if (!histo) {
	MSG_ERROR("Failed to scale histo=NULL in analysis "
		  << name() << " (scale=" << scale << ")");
	return;
      }
      const string hpath = tree().findPath(dynamic_cast<const AIDA::IManagedObject&>(*histo));
      MSG_TRACE("Scaling histo " << hpath);
      
      vector<double> x, y, ex, ey;
      for (size_t i = 0, N = histo->size(); i < N; ++i) {

	IDataPoint * point = histo->point(i);
	assert(point->dimension()==2);
	x .push_back(point->coordinate(0)->value());
	ex.push_back(0.5*(point->coordinate(0)->errorPlus()+
			  point->coordinate(0)->errorMinus()));
	y .push_back(point->coordinate(1)->value()*scale);
	ey.push_back(0.5*scale*(point->coordinate(1)->errorPlus()+
				point->coordinate(1)->errorMinus()));
      }
      string title = histo->title();
      string xtitle = histo->xtitle();
      string ytitle = histo->ytitle();
      
      tree().mkdir("/tmpnormalize");
      tree().mv(hpath, "/tmpnormalize");
      
      if (hpath.find(" ") != string::npos) {
	throw Error("Histogram path '" + hpath + "' is invalid: spaces are not permitted in paths");
      }
      AIDA::IDataPointSet* dps = datapointsetFactory().createXY(hpath, title, x, y, ex, ey);
      dps->setXTitle(xtitle);
      dps->setYTitle(ytitle);
      
      tree().rm(tree().findPath(dynamic_cast<AIDA::IManagedObject&>(*histo)));
      tree().rmdir("/tmpnormalize");
      
      // Set histo pointer to null - it can no longer be used.
      histo = 0;
    }
  };  
      
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SLD_2004_S5693039);

}
