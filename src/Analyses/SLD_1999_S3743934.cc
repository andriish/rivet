// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include "Rivet/Projections/Thrust.hh"
#include "LWH/AIManagedObject.h"
using namespace AIDA;

namespace Rivet {


  /// @brief SLD flavour-dependent fragmentation paper
  /// @author Peter Richardson
  class SLD_1999_S3743934 : public Analysis {
  public:

    /// Constructor
    SLD_1999_S3743934() : Analysis("SLD_1999_S3743934"),
                          _SumOfudsWeights(0.), _SumOfcWeights(0.),
                          _SumOfbWeights(0.),
                          _multPiPlus(4,0.),_multKPlus(4,0.),_multK0(4,0.),
                          _multKStar0(4,0.),_multPhi(4,0.),
                          _multProton(4,0.),_multLambda(4,0.)
    {}


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
      // thrust axis for projections
      Vector3 axis = applyProjection<Thrust>(e, "Thrust").thrustAxis();
      double dot(0.);
      if(!quarks.empty()) {
        dot = quarks[0].momentum().vector3().dot(axis);
        if(quarks[0].pdgId()<0) dot *= -1.;
      }

      foreach (const Particle& p, fs.particles()) {
        const double xp = p.momentum().vector3().mod()/meanBeamMom;
        // if in quark or antiquark hemisphere
        bool quark = p.momentum().vector3().dot(axis)*dot>0.;
        _histXpChargedN->fill(xp, weight);
        int id = abs(p.pdgId());
        // charged pions
        if(id==PIPLUS) {
          _histXpPiPlusN->fill(xp, weight);
          _multPiPlus[0] += weight;
          switch (flavour) {
          case DQUARK: case UQUARK: case SQUARK:
            _multPiPlus[1] += weight;
            _histXpPiPlusLight->fill(xp, weight);
            if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRPiPlus->fill(xp, weight);
            else
              _histRPiMinus->fill(xp, weight);
            break;
          case CQUARK:
            _multPiPlus[2] += weight;
            _histXpPiPlusCharm->fill(xp, weight);
            break;
          case BQUARK:
            _multPiPlus[3] += weight;
            _histXpPiPlusBottom->fill(xp, weight);
            break;
          }
        }
        else if(id==KPLUS) {
          _histXpKPlusN->fill(xp, weight);
          _multKPlus[0] += weight;
          switch (flavour) {
          case DQUARK: case UQUARK: case SQUARK:
            _multKPlus[1] += weight;
            _tempXpKPlusLight->fill(xp, weight);
            _histXpKPlusLight->fill(xp, weight);
            if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRKPlus->fill(xp, weight);
            else
              _histRKMinus->fill(xp, weight);
            break;
         break;
          case CQUARK:
            _multKPlus[2] += weight;
            _histXpKPlusCharm->fill(xp, weight);
            _tempXpKPlusCharm->fill(xp, weight);
            break;
          case BQUARK:
            _multKPlus[3] += weight;
            _histXpKPlusBottom->fill(xp, weight);
            break;
          }
        }
        else if(id==PROTON) {
          _histXpProtonN->fill(xp, weight);
          _multProton[0] += weight;
          switch (flavour) {
          case DQUARK: case UQUARK: case SQUARK:
            _multProton[1] += weight;
            _tempXpProtonLight->fill(xp, weight);
            _histXpProtonLight->fill(xp, weight);
            if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRProton->fill(xp, weight);
            else
              _histRPBar  ->fill(xp, weight);
            break;
         break;
          case CQUARK:
            _multProton[2] += weight;
            _tempXpProtonCharm->fill(xp, weight);
            _histXpProtonCharm->fill(xp, weight);
            break;
          case BQUARK:
            _multProton[3] += weight;
            _histXpProtonBottom->fill(xp, weight);
            break;
          }
        }
      }
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");
      foreach (const Particle& p, ufs.particles()) {
        const double xp = p.momentum().vector3().mod()/meanBeamMom;
        // if in quark or antiquark hemisphere
        bool quark = p.momentum().vector3().dot(axis)*dot>0.;
        int id = abs(p.pdgId());
        if(id==LAMBDA) {
          _multLambda[0] += weight;
          _histXpLambdaN->fill(xp, weight);
          switch (flavour) {
          case DQUARK: case UQUARK: case SQUARK:
            _multLambda[1] += weight;
            _histXpLambdaLight->fill(xp, weight);
            if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRLambda->fill(xp, weight);
            else
              _histRLBar  ->fill(xp, weight);
            break;
          case CQUARK:
            _multLambda[2] += weight;
            _histXpLambdaCharm->fill(xp, weight);
            break;
          case BQUARK:
            _multLambda[3] += weight;
            _histXpLambdaBottom->fill(xp, weight);
            break;
          }
        }
        else if(id==313) {
          _multKStar0[0] += weight;
          _histXpKStar0N->fill(xp, weight);
          switch (flavour) {
          case DQUARK: case UQUARK: case SQUARK:
            _multKStar0[1] += weight;
            _tempXpKStar0Light->fill(xp, weight);
            _histXpKStar0Light->fill(xp, weight);
            if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRKS0   ->fill(xp, weight);
            else
              _histRKSBar0->fill(xp, weight);
            break;
            break;
          case CQUARK:
            _multKStar0[2] += weight;
            _tempXpKStar0Charm->fill(xp, weight);
            _histXpKStar0Charm->fill(xp, weight);
            break;
          case BQUARK:
            _multKStar0[3] += weight;
            _histXpKStar0Bottom->fill(xp, weight);
            break;
          }
        }
        else if(id==333) {
          _multPhi[0] += weight;
          _histXpPhiN->fill(xp, weight);
          switch (flavour) {
          case DQUARK: case UQUARK: case SQUARK:
            _multPhi[1] += weight;
            _histXpPhiLight->fill(xp, weight);
            break;
          case CQUARK:
            _multPhi[2] += weight;
            _histXpPhiCharm->fill(xp, weight);
            break;
          case BQUARK:
            _multPhi[3] += weight;
            _histXpPhiBottom->fill(xp, weight);
            break;
          }
        }
        else if(id==K0S || id == K0L) {
          _multK0[0] += weight;
          _histXpK0N->fill(xp, weight);
          switch (flavour) {
          case DQUARK: case UQUARK: case SQUARK:
            _multK0[1] += weight;
            _histXpK0Light->fill(xp, weight);
            break;
          case CQUARK:
            _multK0[2] += weight;
            _histXpK0Charm->fill(xp, weight);
            break;
          case BQUARK:
            _multK0[3] += weight;
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
      addProjection(Thrust(FinalState()), "Thrust");
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

      _histRPiPlus  = bookHistogram1D( 26, 1, 1);
      _histRPiMinus = bookHistogram1D( 26, 1, 2);
      _histRKS0     = bookHistogram1D( 28, 1, 1);
      _histRKSBar0  = bookHistogram1D( 28, 1, 2);
      _histRKPlus   = bookHistogram1D( 30, 1, 1);
      _histRKMinus  = bookHistogram1D( 30, 1, 2);
      _histRProton  = bookHistogram1D( 32, 1, 1);
      _histRPBar    = bookHistogram1D( 32, 1, 2);
      _histRLambda  = bookHistogram1D( 34, 1, 1);
      _histRLBar    = bookHistogram1D( 34, 1, 2);

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
      // leading particles
      AIDA::IHistogram1D * num = histogramFactory().subtract(dir + "/n1",*_histRPiMinus,*_histRPiPlus);
      AIDA::IHistogram1D * den = histogramFactory().add     (dir + "/n2",*_histRPiMinus,*_histRPiPlus);
      h   = histogramFactory().divide(dir +"/d27-x01-y01",*num,*den);
      histogramFactory().destroy(num);
      histogramFactory().destroy(den);
      num = histogramFactory().subtract(dir + "/n3",*_histRKSBar0,*_histRKS0);
      den = histogramFactory().add     (dir + "/n4",*_histRKSBar0,*_histRKS0);
      h   = histogramFactory().divide(dir +"/d29-x01-y01",*num,*den);
      histogramFactory().destroy(num);
      histogramFactory().destroy(den);
      num = histogramFactory().subtract(dir + "/n5",*_histRKMinus,*_histRKPlus);
      den = histogramFactory().add     (dir + "/n6",*_histRKMinus,*_histRKPlus);
      h   = histogramFactory().divide(dir +"/d31-x01-y01",*num,*den);
      histogramFactory().destroy(num);
      histogramFactory().destroy(den);
      num = histogramFactory().subtract(dir + "/n7",*_histRProton,*_histRPBar);
      den = histogramFactory().add     (dir + "/n8",*_histRProton,*_histRPBar);
      h   = histogramFactory().divide(dir +"/d33-x01-y01",*num,*den);
      histogramFactory().destroy(num);
      histogramFactory().destroy(den);
      num = histogramFactory().subtract(dir + "/n9" ,*_histRLambda,*_histRLBar);
      den = histogramFactory().add     (dir + "/n10",*_histRLambda,*_histRLBar);
      h   = histogramFactory().divide(dir +"/d35-x01-y01",*num,*den);
      histogramFactory().destroy(num);
      histogramFactory().destroy(den);
      // then the rest
      Analysis::scale(_histXpPiPlusN    ,1./sumOfWeights());
      Analysis::scale(_histXpKPlusN     ,1./sumOfWeights());
      Analysis::scale(_histXpProtonN    ,1./sumOfWeights());
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
      Analysis::scale(_histRPiPlus ,1./_SumOfudsWeights);
      Analysis::scale(_histRPiMinus,1./_SumOfudsWeights);
      Analysis::scale(_histRKS0    ,1./_SumOfudsWeights);
      Analysis::scale(_histRKSBar0 ,1./_SumOfudsWeights);
      Analysis::scale(_histRKPlus  ,1./_SumOfudsWeights);
      Analysis::scale(_histRKMinus ,1./_SumOfudsWeights);
      Analysis::scale(_histRProton ,1./_SumOfudsWeights);
      Analysis::scale(_histRPBar   ,1./_SumOfudsWeights);
      Analysis::scale(_histRLambda ,1./_SumOfudsWeights);
      Analysis::scale(_histRLBar   ,1./_SumOfudsWeights);

//       // multiplicities
      AIDA::IDataPointSet * multA;
      AIDA::IDataPointSet * multL;
      AIDA::IDataPointSet * multC;
      AIDA::IDataPointSet * multB;
      AIDA::IDataPointSet * multD1;
      AIDA::IDataPointSet * multD2;
      double  avgNumPartsAll,avgNumPartsLight,avgNumPartsCharm, avgNumPartsBottom;
      // pi+/-
      // all
      avgNumPartsAll = _multPiPlus[0]/sumOfWeights();
      multA = bookDataPointSet(24, 1, 1);
      multA->point(0)->coordinate(1)->setValue(avgNumPartsAll);
      // light
      avgNumPartsLight = _multPiPlus[1]/_SumOfudsWeights;
      multL = bookDataPointSet(24, 1, 2);
      multL->point(0)->coordinate(1)->setValue(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multPiPlus[2]/_SumOfcWeights;
      multC = bookDataPointSet(24, 1, 3);
      multC->point(0)->coordinate(1)->setValue(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multPiPlus[3]/_SumOfbWeights;
      multB = bookDataPointSet(24, 1, 4);
      multB->point(0)->coordinate(1)->setValue(avgNumPartsBottom);
      // charm-light
      multD1 = bookDataPointSet(25, 1, 1);
      multD1->point(0)->coordinate(1)->setValue(avgNumPartsCharm -avgNumPartsLight);
      // bottom-light
      multD2 = bookDataPointSet(25, 1, 2);
      multD2->point(0)->coordinate(1)->setValue(avgNumPartsBottom-avgNumPartsLight);
      // K+/-
      // all
      avgNumPartsAll = _multKPlus[0]/sumOfWeights();
      multA = bookDataPointSet(24, 2, 1);
      multA->point(0)->coordinate(1)->setValue(avgNumPartsAll);
      // light
      avgNumPartsLight = _multKPlus[1]/_SumOfudsWeights;
      multL = bookDataPointSet(24, 2, 2);
      multL->point(0)->coordinate(1)->setValue(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multKPlus[2]/_SumOfcWeights;
      multC = bookDataPointSet(24, 2, 3);
      multC->point(0)->coordinate(1)->setValue(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multKPlus[3]/_SumOfbWeights;
      multB = bookDataPointSet(24, 2, 4);
      multB->point(0)->coordinate(1)->setValue(avgNumPartsBottom);
      // charm-light
      multD1 = bookDataPointSet(25, 2, 1);
      multD1->point(0)->coordinate(1)->setValue(avgNumPartsCharm -avgNumPartsLight);
      // bottom-light
      multD2 = bookDataPointSet(25, 2, 2);
      multD2->point(0)->coordinate(1)->setValue(avgNumPartsBottom-avgNumPartsLight);
      // K0
      // all
      avgNumPartsAll = _multK0[0]/sumOfWeights();
      multA = bookDataPointSet(24, 3, 1);
      multA->point(0)->coordinate(1)->setValue(avgNumPartsAll);
      // light
      avgNumPartsLight = _multK0[1]/_SumOfudsWeights;
      multL = bookDataPointSet(24, 3, 2);
      multL->point(0)->coordinate(1)->setValue(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multK0[2]/_SumOfcWeights;
      multC = bookDataPointSet(24, 3, 3);
      multC->point(0)->coordinate(1)->setValue(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multK0[3]/_SumOfbWeights;
      multB = bookDataPointSet(24, 3, 4);
      multB->point(0)->coordinate(1)->setValue(avgNumPartsBottom);
      // charm-light
      multD1 = bookDataPointSet(25, 3, 1);
      multD1->point(0)->coordinate(1)->setValue(avgNumPartsCharm -avgNumPartsLight);
      // bottom-light
      multD2 = bookDataPointSet(25, 3, 2);
      multD2->point(0)->coordinate(1)->setValue(avgNumPartsBottom-avgNumPartsLight);
      // K*0
      // all
      avgNumPartsAll = _multKStar0[0]/sumOfWeights();
      multA = bookDataPointSet(24, 4, 1);
      multA->point(0)->coordinate(1)->setValue(avgNumPartsAll);
      // light
      avgNumPartsLight = _multKStar0[1]/_SumOfudsWeights;
      multL = bookDataPointSet(24, 4, 2);
      multL->point(0)->coordinate(1)->setValue(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multKStar0[2]/_SumOfcWeights;
      multC = bookDataPointSet(24, 4, 3);
      multC->point(0)->coordinate(1)->setValue(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multKStar0[3]/_SumOfbWeights;
      multB = bookDataPointSet(24, 4, 4);
      multB->point(0)->coordinate(1)->setValue(avgNumPartsBottom);
      // charm-light
      multD1 = bookDataPointSet(25, 4, 1);
      multD1->point(0)->coordinate(1)->setValue(avgNumPartsCharm -avgNumPartsLight);
      // bottom-light
      multD2 = bookDataPointSet(25, 4, 2);
      multD2->point(0)->coordinate(1)->setValue(avgNumPartsBottom-avgNumPartsLight);
      // phi
      // all
      avgNumPartsAll = _multPhi[0]/sumOfWeights();
      multA = bookDataPointSet(24, 5, 1);
      multA->point(0)->coordinate(1)->setValue(avgNumPartsAll);
      // light
      avgNumPartsLight = _multPhi[1]/_SumOfudsWeights;
      multL = bookDataPointSet(24, 5, 2);
      multL->point(0)->coordinate(1)->setValue(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multPhi[2]/_SumOfcWeights;
      multC = bookDataPointSet(24, 5, 3);
      multC->point(0)->coordinate(1)->setValue(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multPhi[3]/_SumOfbWeights;
      multB = bookDataPointSet(24, 5, 4);
      multB->point(0)->coordinate(1)->setValue(avgNumPartsBottom);
      // charm-light
      multD1 = bookDataPointSet(25, 5, 1);
      multD1->point(0)->coordinate(1)->setValue(avgNumPartsCharm -avgNumPartsLight);
      // bottom-light
      multD2 = bookDataPointSet(25, 5, 2);
      multD2->point(0)->coordinate(1)->setValue(avgNumPartsBottom-avgNumPartsLight);
      // p
      // all
      avgNumPartsAll = _multProton[0]/sumOfWeights();
      multA = bookDataPointSet(24, 6, 1);
      multA->point(0)->coordinate(1)->setValue(avgNumPartsAll);
      // light
      avgNumPartsLight = _multProton[1]/_SumOfudsWeights;
      multL = bookDataPointSet(24, 6, 2);
      multL->point(0)->coordinate(1)->setValue(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multProton[2]/_SumOfcWeights;
      multC = bookDataPointSet(24, 6, 3);
      multC->point(0)->coordinate(1)->setValue(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multProton[3]/_SumOfbWeights;
      multB = bookDataPointSet(24, 6, 4);
      multB->point(0)->coordinate(1)->setValue(avgNumPartsBottom);
      // charm-light
      multD1 = bookDataPointSet(25, 6, 1);
      multD1->point(0)->coordinate(1)->setValue(avgNumPartsCharm -avgNumPartsLight);
      // bottom-light
      multD2 = bookDataPointSet(25, 6, 2);
      multD2->point(0)->coordinate(1)->setValue(avgNumPartsBottom-avgNumPartsLight);
      // Lambda
      // all
      avgNumPartsAll = _multLambda[0]/sumOfWeights();
      multA = bookDataPointSet(24, 7, 1);
      multA->point(0)->coordinate(1)->setValue(avgNumPartsAll);
      // light
      avgNumPartsLight = _multLambda[1]/_SumOfudsWeights;
      multL = bookDataPointSet(24, 7, 2);
      multL->point(0)->coordinate(1)->setValue(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multLambda[2]/_SumOfcWeights;
      multC = bookDataPointSet(24, 7, 3);
      multC->point(0)->coordinate(1)->setValue(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multLambda[3]/_SumOfbWeights;
      multB = bookDataPointSet(24, 7, 4);
      multB->point(0)->coordinate(1)->setValue(avgNumPartsBottom);
      // charm-light
      multD1 = bookDataPointSet(25, 7, 1);
      multD1->point(0)->coordinate(1)->setValue(avgNumPartsCharm -avgNumPartsLight);
      // bottom-light
      multD2 = bookDataPointSet(25, 7, 2);
      multD2->point(0)->coordinate(1)->setValue(avgNumPartsBottom-avgNumPartsLight);
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
    vector<double> _multPiPlus,_multKPlus,_multK0,_multKStar0,
      _multPhi,_multProton,_multLambda;

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
    AIDA::IHistogram1D *_histRPiPlus ;
    AIDA::IHistogram1D *_histRPiMinus;
    AIDA::IHistogram1D *_histRKS0    ;
    AIDA::IHistogram1D *_histRKSBar0 ;
    AIDA::IHistogram1D *_histRKPlus  ;
    AIDA::IHistogram1D *_histRKMinus ;
    AIDA::IHistogram1D *_histRProton ;
    AIDA::IHistogram1D *_histRPBar   ;
    AIDA::IHistogram1D *_histRLambda ;
    AIDA::IHistogram1D *_histRLBar   ;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SLD_1999_S3743934);

}
