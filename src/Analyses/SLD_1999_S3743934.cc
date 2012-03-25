// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include "LWH/AIManagedObject.h"
using namespace AIDA;

namespace Rivet {


  /// @brief SLD flavour-dependent fragmentation paper
  /// @author Peter Richardson
  class SLD_1999_S3743934 : public Analysis {
  public:

    /// Constructor
    SLD_1999_S3743934() : Analysis("SLD_1999_S3743934")
    {
      // Counters
      _SumOfudsWeights = 0;
      _SumOfcWeights = 0;
      _SumOfbWeights = 0;
    }


    /// @name Analysis methods
    //@{

    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
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
      switch (flavour) {
      case 1: case 2: case 3:
        _SumOfudsWeights += weight;
        break;
      case 4:
        _SumOfcWeights += weight;
        break;
      case 5:
        _SumOfbWeights += weight;
        break;
      }

      foreach (const Particle& p, fs.particles()) {
	const double xp = p.momentum().vector3().mod()/meanBeamMom;
	_histXpChargedN->fill(xp, weight);
	int id = abs(p.pdgId());
	// charged pions
	if(id==PIPLUS) {
	  _histXpPiPlusN->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _histXpPiPlusLight->fill(xp, weight);
         break;
	  case CQUARK:
	    _histXpPiPlusCharm->fill(xp, weight);
	    break;
	  case BQUARK:
	    _histXpPiPlusBottom->fill(xp, weight);
	    break;
	  }
	}
	else if(id==KPLUS) {
	  _histXpKPlusN->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _tempXpKPlusLight->fill(xp, weight);
	    _histXpKPlusLight->fill(xp, weight);
         break;
	  case CQUARK:
	    _histXpKPlusCharm->fill(xp, weight);
	    _tempXpKPlusCharm->fill(xp, weight);
	    break;
	  case BQUARK:
	    _histXpKPlusBottom->fill(xp, weight);
	    break;
	  }
	}
	else if(id==PROTON) {
	  _histXpProtonN->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _tempXpProtonLight->fill(xp, weight);
	    _histXpProtonLight->fill(xp, weight);
         break;
	  case CQUARK:
	    _tempXpProtonCharm->fill(xp, weight);
	    _histXpProtonCharm->fill(xp, weight);
	    break;
	  case BQUARK:
	    _histXpProtonBottom->fill(xp, weight);
	    break;
	  }
	}
      }
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");
      foreach (const Particle& p, ufs.particles()) {
	const double xp = p.momentum().vector3().mod()/meanBeamMom;
	int id = abs(p.pdgId());
	if(id==LAMBDA) {
	  _histXpLambdaN->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _histXpLambdaLight->fill(xp, weight);
	    break;
	  case CQUARK:
	    _histXpLambdaCharm->fill(xp, weight);
	    break;
	  case BQUARK:
	    _histXpLambdaBottom->fill(xp, weight);
	    break;
	  }
	}
	else if(id==313) {
	  _histXpKStar0N->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _tempXpKStar0Light->fill(xp, weight);
	    _histXpKStar0Light->fill(xp, weight);
	    break;
	  case CQUARK:
	    _tempXpKStar0Charm->fill(xp, weight);
	    _histXpKStar0Charm->fill(xp, weight);
	    break;
	  case BQUARK:
	    _histXpKStar0Bottom->fill(xp, weight);
	    break;
	  }
	}
	else if(id==333) {
	  _histXpPhiN->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _histXpPhiLight->fill(xp, weight);
	    break;
	  case CQUARK:
	    _histXpPhiCharm->fill(xp, weight);
	    break;
	  case BQUARK:
	    _histXpPhiBottom->fill(xp, weight);
	    break;
	  }
	}
	else if(id==K0S || id == K0L) {
	  _histXpK0N->fill(xp, weight);
	  switch (flavour) {
	  case DQUARK: case UQUARK: case SQUARK:
	    _histXpK0Light->fill(xp, weight);
	    break;
	  case CQUARK:
	    _histXpK0Charm->fill(xp, weight);
	    break;
	  case BQUARK:
	    _histXpK0Bottom->fill(xp, weight);
	    break;
	  }
	}
      }
    }


    void init() {
      // Projections
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "FS");
      addProjection(UnstableFinalState(), "UFS");
      addProjection(InitialQuarks(), "IQF");

      _histXpPiPlusN      = bookHistogram1D( 1, 1, 2);
      _histXpKPlusN       = bookHistogram1D( 2, 1, 2);
      _histXpProtonN      = bookHistogram1D( 3, 1, 2);
      _histXpChargedN     = bookHistogram1D( 4, 1, 1);
      _histXpK0N          = bookHistogram1D( 5, 1, 1);
      _histXpLambdaN      = bookHistogram1D( 7, 1, 1);
      _histXpKStar0N      = bookHistogram1D( 8, 1, 1);
      _histXpPhiN         = bookHistogram1D( 9, 1, 1);
      _histXpPiPlusLight  = bookHistogram1D(10, 1, 1);
      _histXpPiPlusCharm  = bookHistogram1D(10, 1, 2);
      _histXpPiPlusBottom = bookHistogram1D(10, 1, 3);
      _histXpKPlusLight   = bookHistogram1D(12, 1, 1);
      _histXpKPlusCharm   = bookHistogram1D(12, 1, 2);
      _histXpKPlusBottom  = bookHistogram1D(12, 1, 3);
      _histXpKStar0Light   = bookHistogram1D(14, 1, 1);
      _histXpKStar0Charm   = bookHistogram1D(14, 1, 2);
      _histXpKStar0Bottom  = bookHistogram1D(14, 1, 3);
      _histXpProtonLight   = bookHistogram1D(16, 1, 1);
      _histXpProtonCharm   = bookHistogram1D(16, 1, 2);
      _histXpProtonBottom  = bookHistogram1D(16, 1, 3);
      _histXpLambdaLight   = bookHistogram1D(18, 1, 1);
      _histXpLambdaCharm   = bookHistogram1D(18, 1, 2);
      _histXpLambdaBottom  = bookHistogram1D(18, 1, 3);
      _histXpK0Light   = bookHistogram1D(20, 1, 1);
      _histXpK0Charm   = bookHistogram1D(20, 1, 2);
      _histXpK0Bottom  = bookHistogram1D(20, 1, 3);
      _histXpPhiLight   = bookHistogram1D(22, 1, 1);
      _histXpPhiCharm   = bookHistogram1D(22, 1, 2);
      _histXpPhiBottom  = bookHistogram1D(22, 1, 3);
      _tempXpKPlusCharm  = bookHistogram1D( "tempXpKPlusCharm" ,
					    binEdges(13,1,1),"tempXpKPlusCharm" ,"X","Y");
      _tempXpKPlusLight  = bookHistogram1D( "tempXpKPlusLight" ,
					    binEdges(13,1,1),"tempXpKPlusLight" ,"X","Y");
      _tempXpKStar0Charm = bookHistogram1D( "tempXpKStar0Charm",
					    binEdges(15,1,1),"tempXpKStar0Charm","X","Y");
      _tempXpKStar0Light = bookHistogram1D( "tempXpKStar0Light",
					    binEdges(15,1,1),"tempXpKStar0Light","X","Y");
      _tempXpProtonCharm = bookHistogram1D( "tempXpProtonCharm",
					    binEdges(17,1,1),"tempXpProtonCharm","X","Y");
      _tempXpProtonLight = bookHistogram1D( "tempXpProtonLight",
					    binEdges(17,1,1),"tempXpProtonLight","X","Y");
    }


    /// Finalize
    void finalize() {
      // get the ratio plots sorted out first
      AIDA::IDataPointSet* h = 0;
      const string dir = histoDir();
      h = histogramFactory().divide(dir + "/d01-x01-y01", *_histXpPiPlusN , *_histXpChargedN );
      h = histogramFactory().divide(dir + "/d02-x01-y01", *_histXpKPlusN , *_histXpChargedN );
      h = histogramFactory().divide(dir + "/d03-x01-y01", *_histXpProtonN , *_histXpChargedN );
      h = histogramFactory().divide(dir + "/d11-x01-y01", *_histXpPiPlusCharm , *_histXpPiPlusLight);
      scale(h,_SumOfudsWeights/_SumOfcWeights);
      h = histogramFactory().divide(dir + "/d11-x01-y02", *_histXpPiPlusBottom, *_histXpPiPlusLight);
      scale(h,_SumOfudsWeights/_SumOfbWeights);
      h = histogramFactory().divide(dir + "/d13-x01-y01", *_tempXpKPlusCharm  , *_tempXpKPlusLight );
      scale(h,_SumOfudsWeights/_SumOfcWeights);
      h = histogramFactory().divide(dir + "/d13-x01-y02", *_histXpKPlusBottom , *_histXpKPlusLight );
      scale(h,_SumOfudsWeights/_SumOfbWeights);
      h = histogramFactory().divide(dir + "/d15-x01-y01", *_tempXpKStar0Charm , *_tempXpKStar0Light);
      scale(h,_SumOfudsWeights/_SumOfcWeights);
      h = histogramFactory().divide(dir + "/d15-x01-y02", *_histXpKStar0Bottom, *_histXpKStar0Light);
      scale(h,_SumOfudsWeights/_SumOfbWeights);
      h = histogramFactory().divide(dir + "/d17-x01-y01", *_tempXpProtonCharm , *_tempXpProtonLight);
      scale(h,_SumOfudsWeights/_SumOfcWeights);
      h = histogramFactory().divide(dir + "/d17-x01-y02", *_histXpProtonBottom, *_histXpProtonLight);
      scale(h,_SumOfudsWeights/_SumOfbWeights);
      h = histogramFactory().divide(dir + "/d19-x01-y01", *_histXpLambdaCharm , *_histXpLambdaLight);
      scale(h,_SumOfudsWeights/_SumOfcWeights);
      h = histogramFactory().divide(dir + "/d19-x01-y02", *_histXpLambdaBottom, *_histXpLambdaLight);
      scale(h,_SumOfudsWeights/_SumOfbWeights);
      h = histogramFactory().divide(dir + "/d21-x01-y01", *_histXpK0Charm     , *_histXpK0Light    );
      scale(h,_SumOfudsWeights/_SumOfcWeights);
      h = histogramFactory().divide(dir + "/d21-x01-y02", *_histXpK0Bottom    , *_histXpK0Light    );
      scale(h,_SumOfudsWeights/_SumOfbWeights);
      h = histogramFactory().divide(dir + "/d23-x01-y01", *_histXpPhiCharm    , *_histXpPhiLight   );
      scale(h,_SumOfudsWeights/_SumOfcWeights);
      h =  histogramFactory().divide(dir + "/d23-x01-y02", *_histXpPhiBottom   , *_histXpPhiLight   );
      scale(h,_SumOfudsWeights/_SumOfbWeights);
      histogramFactory().destroy(_tempXpKPlusCharm );
      histogramFactory().destroy(_tempXpKPlusLight );
      histogramFactory().destroy(_tempXpKStar0Charm);
      histogramFactory().destroy(_tempXpKStar0Light);
      histogramFactory().destroy(_tempXpProtonCharm);
      histogramFactory().destroy(_tempXpProtonLight);
      // then the rest
      Analysis::scale(_histXpPiPlusN,1./sumOfWeights());
      Analysis::scale(_histXpKPlusN,1./sumOfWeights());
      Analysis::scale(_histXpProtonN,1./sumOfWeights());
      Analysis::scale(_histXpChargedN,1./sumOfWeights());
      Analysis::scale(_histXpK0N,1./sumOfWeights());
      Analysis::scale(_histXpLambdaN,1./sumOfWeights());
      Analysis::scale(_histXpKStar0N,1./sumOfWeights());
      Analysis::scale(_histXpPhiN,1./sumOfWeights());
      Analysis::scale(_histXpPiPlusLight,1./_SumOfudsWeights);
      Analysis::scale(_histXpPiPlusCharm,1./_SumOfcWeights);
      Analysis::scale(_histXpPiPlusBottom,1./_SumOfbWeights);
      Analysis::scale(_histXpKPlusLight ,1./_SumOfudsWeights);
      Analysis::scale(_histXpKPlusCharm ,1./_SumOfcWeights);
      Analysis::scale(_histXpKPlusBottom,1./_SumOfbWeights);
      Analysis::scale(_histXpKStar0Light ,1./_SumOfudsWeights);
      Analysis::scale(_histXpKStar0Charm ,1./_SumOfcWeights);
      Analysis::scale(_histXpKStar0Bottom,1./_SumOfbWeights);
      Analysis::scale(_histXpProtonLight ,1./_SumOfudsWeights);
      Analysis::scale(_histXpProtonCharm ,1./_SumOfcWeights);
      Analysis::scale(_histXpProtonBottom,1./_SumOfbWeights);
      Analysis::scale(_histXpLambdaLight ,1./_SumOfudsWeights);
      Analysis::scale(_histXpLambdaCharm ,1./_SumOfcWeights);
      Analysis::scale(_histXpLambdaBottom,1./_SumOfbWeights);
      Analysis::scale(_histXpK0Light ,1./_SumOfudsWeights);
      Analysis::scale(_histXpK0Charm ,1./_SumOfcWeights);
      Analysis::scale(_histXpK0Bottom,1./_SumOfbWeights);
      Analysis::scale(_histXpPhiLight ,1./_SumOfudsWeights);
      Analysis::scale(_histXpPhiCharm ,1./_SumOfcWeights);
      Analysis::scale(_histXpPhiBottom,1./_SumOfbWeights);
    }

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

  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.
    double _SumOfudsWeights,_SumOfcWeights,_SumOfbWeights;

    AIDA::IHistogram1D *_histXpPiPlusSig;   
    AIDA::IHistogram1D *_histXpPiPlusN;     
    AIDA::IHistogram1D *_histXpKPlusSig;    
    AIDA::IHistogram1D *_histXpKPlusN;      
    AIDA::IHistogram1D *_histXpProtonSig;   
    AIDA::IHistogram1D *_histXpProtonN;     
    AIDA::IHistogram1D *_histXpChargedN;    
    AIDA::IHistogram1D *_histXpK0N;    
    AIDA::IHistogram1D *_histXpLambdaN;     
    AIDA::IHistogram1D *_histXpKStar0N;     
    AIDA::IHistogram1D *_histXpPhiN;       
    AIDA::IHistogram1D *_histXpPiPlusLight; 
    AIDA::IHistogram1D *_histXpPiPlusCharm; 
    AIDA::IHistogram1D *_histXpPiPlusBottom;
    AIDA::IHistogram1D *_histXpKPlusLight;  
    AIDA::IHistogram1D *_histXpKPlusCharm;  
    AIDA::IHistogram1D *_histXpKPlusBottom;
    AIDA::IHistogram1D *_histXpKStar0Light;   
    AIDA::IHistogram1D *_histXpKStar0Charm;   
    AIDA::IHistogram1D *_histXpKStar0Bottom;
    AIDA::IHistogram1D *_histXpProtonLight;   
    AIDA::IHistogram1D *_histXpProtonCharm;   
    AIDA::IHistogram1D *_histXpProtonBottom;
    AIDA::IHistogram1D *_histXpLambdaLight;   
    AIDA::IHistogram1D *_histXpLambdaCharm;   
    AIDA::IHistogram1D *_histXpLambdaBottom;
    AIDA::IHistogram1D *_histXpK0Light;
    AIDA::IHistogram1D *_histXpK0Charm;
    AIDA::IHistogram1D *_histXpK0Bottom;
    AIDA::IHistogram1D *_histXpPhiLight;
    AIDA::IHistogram1D *_histXpPhiCharm;
    AIDA::IHistogram1D *_histXpPhiBottom;
    AIDA::IHistogram1D *_tempXpKPlusCharm ;
    AIDA::IHistogram1D *_tempXpKPlusLight ;
    AIDA::IHistogram1D *_tempXpKStar0Charm;
    AIDA::IHistogram1D *_tempXpKStar0Light;
    AIDA::IHistogram1D *_tempXpProtonCharm;
    AIDA::IHistogram1D *_tempXpProtonLight;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SLD_1999_S3743934);

}
