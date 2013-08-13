// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief SLD flavour-dependent fragmentation paper
  /// @author Peter Richardson
  class SLD_1999_S3743934 : public Analysis {
  public:

    /// Constructor
    SLD_1999_S3743934()
      : Analysis("SLD_1999_S3743934"),
        _SumOfudsWeights(0.), _SumOfcWeights(0.),
        _SumOfbWeights(0.),
        _multPiPlus(4,0.),_multKPlus(4,0.),_multK0(4,0.),
        _multKStar0(4,0.),_multPhi(4,0.),
        _multProton(4,0.),_multLambda(4,0.)
    {    }


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
      /// @todo Can we make this based on hadron flavour instead?
      Particles quarks;
      if (iqf.particles().size() == 2) {
        flavour = abs( iqf.particles().front().pdgId() );
        quarks = iqf.particles();
      } else {
        map<int, Particle > quarkmap;
        foreach (const Particle& p, iqf.particles()) {
          if (quarkmap.find(p.pdgId()) == quarkmap.end()) quarkmap[p.pdgId()] = p;
          else if (quarkmap[p.pdgId()].momentum().E() < p.momentum().E()) quarkmap[p.pdgId()] = p;
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          double energy(0.);
          if (quarkmap.find( i) != quarkmap.end())
            energy += quarkmap[ i].momentum().E();
          if (quarkmap.find(-i) != quarkmap.end())
            energy += quarkmap[-i].momentum().E();
          if (energy > maxenergy)
            flavour = i;
        }
        if (quarkmap.find(flavour) != quarkmap.end())
          quarks.push_back(quarkmap[flavour]);
        if (quarkmap.find(-flavour) != quarkmap.end())
          quarks.push_back(quarkmap[-flavour]);
      }
      switch (flavour) {
      case PID::DQUARK:
      case PID::UQUARK:
      case PID::SQUARK:
        _SumOfudsWeights += weight;
        break;
      case PID::CQUARK:
        _SumOfcWeights += weight;
        break;
      case PID::BQUARK:
        _SumOfbWeights += weight;
        break;
      }
      // thrust axis for projections
      Vector3 axis = applyProjection<Thrust>(e, "Thrust").thrustAxis();
      double dot(0.);
      if (!quarks.empty()) {
        dot = quarks[0].momentum().vector3().dot(axis);
        if (quarks[0].pdgId() < 0) dot *= -1;
      }

      foreach (const Particle& p, fs.particles()) {
        const double xp = p.momentum().vector3().mod()/meanBeamMom;
        // if in quark or antiquark hemisphere
        bool quark = p.momentum().vector3().dot(axis)*dot>0.;
        _histXpChargedN->fill(xp, weight);
        int id = abs(p.pdgId());
        // charged pions
        if (id == PID::PIPLUS) {
          _histXpPiPlusN->fill(xp, weight);
          _multPiPlus[0] += weight;
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multPiPlus[1] += weight;
            _histXpPiPlusLight->fill(xp, weight);
            if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRPiPlus->fill(xp, weight);
            else
              _histRPiMinus->fill(xp, weight);
            break;
          case PID::CQUARK:
            _multPiPlus[2] += weight;
            _histXpPiPlusCharm->fill(xp, weight);
            break;
          case PID::BQUARK:
            _multPiPlus[3] += weight;
            _histXpPiPlusBottom->fill(xp, weight);
            break;
          }
        }
        else if (id == PID::KPLUS) {
          _histXpKPlusN->fill(xp, weight);
          _multKPlus[0] += weight;
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multKPlus[1] += weight;
            _tempXpKPlusLight->fill(xp, weight);
            _histXpKPlusLight->fill(xp, weight);
            if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRKPlus->fill(xp, weight);
            else
              _histRKMinus->fill(xp, weight);
            break;
         break;
          case PID::CQUARK:
            _multKPlus[2] += weight;
            _histXpKPlusCharm->fill(xp, weight);
            _tempXpKPlusCharm->fill(xp, weight);
            break;
          case PID::BQUARK:
            _multKPlus[3] += weight;
            _histXpKPlusBottom->fill(xp, weight);
            break;
          }
        }
        else if (id == PID::PROTON) {
          _histXpProtonN->fill(xp, weight);
          _multProton[0] += weight;
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multProton[1] += weight;
            _tempXpProtonLight->fill(xp, weight);
            _histXpProtonLight->fill(xp, weight);
            if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRProton->fill(xp, weight);
            else
              _histRPBar  ->fill(xp, weight);
            break;
         break;
          case PID::CQUARK:
            _multProton[2] += weight;
            _tempXpProtonCharm->fill(xp, weight);
            _histXpProtonCharm->fill(xp, weight);
            break;
          case PID::BQUARK:
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
        if (id == PID::LAMBDA) {
          _multLambda[0] += weight;
          _histXpLambdaN->fill(xp, weight);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multLambda[1] += weight;
            _histXpLambdaLight->fill(xp, weight);
            if( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRLambda->fill(xp, weight);
            else
              _histRLBar  ->fill(xp, weight);
            break;
          case PID::CQUARK:
            _multLambda[2] += weight;
            _histXpLambdaCharm->fill(xp, weight);
            break;
          case PID::BQUARK:
            _multLambda[3] += weight;
            _histXpLambdaBottom->fill(xp, weight);
            break;
          }
        }
        else if (id == 313) {
          _multKStar0[0] += weight;
          _histXpKStar0N->fill(xp, weight);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multKStar0[1] += weight;
            _tempXpKStar0Light->fill(xp, weight);
            _histXpKStar0Light->fill(xp, weight);
            if ( ( quark && p.pdgId()>0 ) || ( !quark && p.pdgId()<0 ))
              _histRKS0   ->fill(xp, weight);
            else
              _histRKSBar0->fill(xp, weight);
            break;
            break;
          case PID::CQUARK:
            _multKStar0[2] += weight;
            _tempXpKStar0Charm->fill(xp, weight);
            _histXpKStar0Charm->fill(xp, weight);
            break;
          case PID::BQUARK:
            _multKStar0[3] += weight;
            _histXpKStar0Bottom->fill(xp, weight);
            break;
          }
        }
        else if (id == 333) {
          _multPhi[0] += weight;
          _histXpPhiN->fill(xp, weight);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multPhi[1] += weight;
            _histXpPhiLight->fill(xp, weight);
            break;
          case PID::CQUARK:
            _multPhi[2] += weight;
            _histXpPhiCharm->fill(xp, weight);
            break;
          case PID::BQUARK:
            _multPhi[3] += weight;
            _histXpPhiBottom->fill(xp, weight);
            break;
          }
        }
        else if (id == PID::K0S || id == PID::K0L) {
          _multK0[0] += weight;
          _histXpK0N->fill(xp, weight);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multK0[1] += weight;
            _histXpK0Light->fill(xp, weight);
            break;
          case PID::CQUARK:
            _multK0[2] += weight;
            _histXpK0Charm->fill(xp, weight);
            break;
          case PID::BQUARK:
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

      _histXpPiPlusN      = bookHisto1D( 1, 1, 2);
      _histXpKPlusN       = bookHisto1D( 2, 1, 2);
      _histXpProtonN      = bookHisto1D( 3, 1, 2);
      _histXpChargedN     = bookHisto1D( 4, 1, 1);
      _histXpK0N          = bookHisto1D( 5, 1, 1);
      _histXpLambdaN      = bookHisto1D( 7, 1, 1);
      _histXpKStar0N      = bookHisto1D( 8, 1, 1);
      _histXpPhiN         = bookHisto1D( 9, 1, 1);

      _histXpPiPlusLight  = bookHisto1D(10, 1, 1);
      _histXpPiPlusCharm  = bookHisto1D(10, 1, 2);
      _histXpPiPlusBottom = bookHisto1D(10, 1, 3);
      _histXpKPlusLight   = bookHisto1D(12, 1, 1);
      _histXpKPlusCharm   = bookHisto1D(12, 1, 2);
      _histXpKPlusBottom  = bookHisto1D(12, 1, 3);
      _histXpKStar0Light  = bookHisto1D(14, 1, 1);
      _histXpKStar0Charm  = bookHisto1D(14, 1, 2);
      _histXpKStar0Bottom = bookHisto1D(14, 1, 3);
      _histXpProtonLight  = bookHisto1D(16, 1, 1);
      _histXpProtonCharm  = bookHisto1D(16, 1, 2);
      _histXpProtonBottom = bookHisto1D(16, 1, 3);
      _histXpLambdaLight  = bookHisto1D(18, 1, 1);
      _histXpLambdaCharm  = bookHisto1D(18, 1, 2);
      _histXpLambdaBottom = bookHisto1D(18, 1, 3);
      _histXpK0Light      = bookHisto1D(20, 1, 1);
      _histXpK0Charm      = bookHisto1D(20, 1, 2);
      _histXpK0Bottom     = bookHisto1D(20, 1, 3);
      _histXpPhiLight     = bookHisto1D(22, 1, 1);
      _histXpPhiCharm     = bookHisto1D(22, 1, 2);
      _histXpPhiBottom    = bookHisto1D(22, 1, 3);

      _tempXpKPlusCharm   = bookHisto1D(13, 1, 1, "tempXpKPlusCharm");
      _tempXpKPlusLight   = bookHisto1D(13, 1, 1, "tempXpKPlusLight");
      _tempXpKStar0Charm  = bookHisto1D(15, 1, 1, "tempXpKStar0Charm");
      _tempXpKStar0Light  = bookHisto1D(15, 1, 1, "tempXpKStar0Light");
      _tempXpProtonCharm  = bookHisto1D(17, 1, 1, "tempXpProtonCharm");
      _tempXpProtonLight  = bookHisto1D(17, 1, 1, "tempXpProtonLight");

      _histRPiPlus  = bookHisto1D( 26, 1, 1);
      _histRPiMinus = bookHisto1D( 26, 1, 2);
      _histRKS0     = bookHisto1D( 28, 1, 1);
      _histRKSBar0  = bookHisto1D( 28, 1, 2);
      _histRKPlus   = bookHisto1D( 30, 1, 1);
      _histRKMinus  = bookHisto1D( 30, 1, 2);
      _histRProton  = bookHisto1D( 32, 1, 1);
      _histRPBar    = bookHisto1D( 32, 1, 2);
      _histRLambda  = bookHisto1D( 34, 1, 1);
      _histRLBar    = bookHisto1D( 34, 1, 2);

      _h_Xp_PiPl_Ch        = bookScatter2D(1, 1, 1);
      _h_Xp_KPl_Ch         = bookScatter2D(2, 1, 1);
      _h_Xp_Pr_Ch          = bookScatter2D(3, 1, 1);
      _h_Xp_PiPlCh_PiPlLi  = bookScatter2D(11, 1, 1);
      _h_Xp_PiPlBo_PiPlLi  = bookScatter2D(11, 1, 2);
      _h_Xp_KPlCh_KPlLi    = bookScatter2D(13, 1, 1);
      _h_Xp_KPlBo_KPlLi    = bookScatter2D(13, 1, 2);
      _h_Xp_KS0Ch_KS0Li    = bookScatter2D(15, 1, 1);
      _h_Xp_KS0Bo_KS0Li    = bookScatter2D(15, 1, 2);
      _h_Xp_PrCh_PrLi      = bookScatter2D(17, 1, 1);
      _h_Xp_PrBo_PrLi      = bookScatter2D(17, 1, 2);
      _h_Xp_LaCh_LaLi      = bookScatter2D(19, 1, 1);
      _h_Xp_LaBo_LaLi      = bookScatter2D(19, 1, 2);
      _h_Xp_K0Ch_K0Li      = bookScatter2D(21, 1, 1);
      _h_Xp_K0Bo_K0Li      = bookScatter2D(21, 1, 2);
      _h_Xp_PhiCh_PhiLi    = bookScatter2D(23, 1, 1);
      _h_Xp_PhiBo_PhiLi    = bookScatter2D(23, 1, 2);

      _h_PiM_PiP    = bookScatter2D(27, 1, 1);
      _h_KSBar0_KS0 = bookScatter2D(29, 1, 1);
      _h_KM_KP      = bookScatter2D(31, 1, 1);
      _h_Pr_PBar    = bookScatter2D(33, 1, 1);
      _h_Lam_LBar   = bookScatter2D(35, 1, 1);
    }


    /// Finalize
    void finalize() {
      // Get the ratio plots sorted out first
      divide(_histXpPiPlusN, _histXpChargedN, _h_Xp_PiPl_Ch);
      divide(_histXpKPlusN, _histXpChargedN, _h_Xp_KPl_Ch);
      divide(_histXpProtonN, _histXpChargedN, _h_Xp_Pr_Ch);
      divide(_histXpPiPlusCharm, _histXpPiPlusLight, _h_Xp_PiPlCh_PiPlLi);
      divide(_histXpPiPlusBottom, _histXpPiPlusLight, _h_Xp_PiPlBo_PiPlLi);
      divide(_tempXpKPlusCharm , _tempXpKPlusLight, _h_Xp_KPlCh_KPlLi);
      divide(_histXpKPlusBottom, _histXpKPlusLight, _h_Xp_KPlBo_KPlLi);
      divide(_tempXpKStar0Charm, _tempXpKStar0Light, _h_Xp_KS0Ch_KS0Li);
      divide(_histXpKStar0Bottom, _histXpKStar0Light, _h_Xp_KS0Bo_KS0Li);
      divide(_tempXpProtonCharm, _tempXpProtonLight, _h_Xp_PrCh_PrLi);
      divide(_histXpProtonBottom, _histXpProtonLight, _h_Xp_PrBo_PrLi);
      divide(_histXpLambdaCharm, _histXpLambdaLight, _h_Xp_LaCh_LaLi);
      divide(_histXpLambdaBottom, _histXpLambdaLight, _h_Xp_LaBo_LaLi);
      divide(_histXpK0Charm, _histXpK0Light, _h_Xp_K0Ch_K0Li);
      divide(_histXpK0Bottom, _histXpK0Light, _h_Xp_K0Bo_K0Li);
      divide(_histXpPhiCharm, _histXpPhiLight, _h_Xp_PhiCh_PhiLi);
      divide(_histXpPhiBottom, _histXpPhiLight, _h_Xp_PhiBo_PhiLi);

      // Then the leading particles
      // divide(*_histRPiMinus - *_histRPiPlus, *_histRPiMinus + *_histRPiPlus, _h_PiM_PiP);
      // divide(*_histRKSBar0 - *_histRKS0, *_histRKSBar0 + *_histRKS0, _h_KSBar0_KS0);
      // divide(*_histRKMinus - *_histRKPlus, *_histRKMinus + *_histRKPlus, _h_KM_KP);
      // divide(*_histRProton - *_histRPBar, *_histRProton + *_histRPBar, _h_Pr_PBar);
      // divide(*_histRLambda - *_histRLBar, *_histRLambda + *_histRLBar, _h_Lam_LBar);

      // Then the rest
      scale(_histXpPiPlusN,      1/sumOfWeights());
      scale(_histXpKPlusN,       1/sumOfWeights());
      scale(_histXpProtonN,      1/sumOfWeights());
      scale(_histXpChargedN,     1/sumOfWeights());
      scale(_histXpK0N,          1/sumOfWeights());
      scale(_histXpLambdaN,      1/sumOfWeights());
      scale(_histXpKStar0N,      1/sumOfWeights());
      scale(_histXpPhiN,         1/sumOfWeights());
      scale(_histXpPiPlusLight,  1/_SumOfudsWeights);
      scale(_histXpPiPlusCharm,  1/_SumOfcWeights);
      scale(_histXpPiPlusBottom, 1/_SumOfbWeights);
      scale(_histXpKPlusLight,   1/_SumOfudsWeights);
      scale(_histXpKPlusCharm,   1/_SumOfcWeights);
      scale(_histXpKPlusBottom,  1/_SumOfbWeights);
      scale(_histXpKStar0Light,  1/_SumOfudsWeights);
      scale(_histXpKStar0Charm,  1/_SumOfcWeights);
      scale(_histXpKStar0Bottom, 1/_SumOfbWeights);
      scale(_histXpProtonLight,  1/_SumOfudsWeights);
      scale(_histXpProtonCharm,  1/_SumOfcWeights);
      scale(_histXpProtonBottom, 1/_SumOfbWeights);
      scale(_histXpLambdaLight,  1/_SumOfudsWeights);
      scale(_histXpLambdaCharm,  1/_SumOfcWeights);
      scale(_histXpLambdaBottom, 1/_SumOfbWeights);
      scale(_histXpK0Light,      1/_SumOfudsWeights);
      scale(_histXpK0Charm,      1/_SumOfcWeights);
      scale(_histXpK0Bottom,     1/_SumOfbWeights);
      scale(_histXpPhiLight,     1/_SumOfudsWeights);
      scale(_histXpPhiCharm ,    1/_SumOfcWeights);
      scale(_histXpPhiBottom,    1/_SumOfbWeights);
      scale(_histRPiPlus,        1/_SumOfudsWeights);
      scale(_histRPiMinus,       1/_SumOfudsWeights);
      scale(_histRKS0,           1/_SumOfudsWeights);
      scale(_histRKSBar0,        1/_SumOfudsWeights);
      scale(_histRKPlus,         1/_SumOfudsWeights);
      scale(_histRKMinus,        1/_SumOfudsWeights);
      scale(_histRProton,        1/_SumOfudsWeights);
      scale(_histRPBar,          1/_SumOfudsWeights);
      scale(_histRLambda,        1/_SumOfudsWeights);
      scale(_histRLBar,          1/_SumOfudsWeights);

      // Multiplicities
      double avgNumPartsAll, avgNumPartsLight,avgNumPartsCharm, avgNumPartsBottom;
      // pi+/-
      // all
      avgNumPartsAll = _multPiPlus[0]/sumOfWeights();
      bookScatter2D(24, 1, 1, true)->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = _multPiPlus[1]/_SumOfudsWeights;
      bookScatter2D(24, 1, 2, true)->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multPiPlus[2]/_SumOfcWeights;
      bookScatter2D(24, 1, 3, true)->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multPiPlus[3]/_SumOfbWeights;
      bookScatter2D(24, 1, 4, true)->point(0).setY(avgNumPartsBottom);
      // charm-light
      bookScatter2D(25, 1, 1, true)->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      bookScatter2D(25, 1, 2, true)->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // K+/-
      // all
      avgNumPartsAll = _multKPlus[0]/sumOfWeights();
      bookScatter2D(24, 2, 1, true)->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = _multKPlus[1]/_SumOfudsWeights;
      bookScatter2D(24, 2, 2, true)->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multKPlus[2]/_SumOfcWeights;
      bookScatter2D(24, 2, 3, true)->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multKPlus[3]/_SumOfbWeights;
      bookScatter2D(24, 2, 4, true)->point(0).setY(avgNumPartsBottom);
      // charm-light
      bookScatter2D(25, 2, 1, true)->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      bookScatter2D(25, 2, 2, true)->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // K0
      // all
      avgNumPartsAll = _multK0[0]/sumOfWeights();
      bookScatter2D(24, 3, 1, true)->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = _multK0[1]/_SumOfudsWeights;
      bookScatter2D(24, 3, 2, true)->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multK0[2]/_SumOfcWeights;
      bookScatter2D(24, 3, 3, true)->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multK0[3]/_SumOfbWeights;
      bookScatter2D(24, 3, 4, true)->point(0).setY(avgNumPartsBottom);
      // charm-light
      bookScatter2D(25, 3, 1, true)->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      bookScatter2D(25, 3, 2, true)->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // K*0
      // all
      avgNumPartsAll = _multKStar0[0]/sumOfWeights();
      bookScatter2D(24, 4, 1, true)->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = _multKStar0[1]/_SumOfudsWeights;
      bookScatter2D(24, 4, 2, true)->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multKStar0[2]/_SumOfcWeights;
      bookScatter2D(24, 4, 3, true)->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multKStar0[3]/_SumOfbWeights;
      bookScatter2D(24, 4, 4, true)->point(0).setY(avgNumPartsBottom);
      // charm-light
      bookScatter2D(25, 4, 1, true)->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      bookScatter2D(25, 4, 2, true)->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // phi
      // all
      avgNumPartsAll = _multPhi[0]/sumOfWeights();
      bookScatter2D(24, 5, 1, true)->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = _multPhi[1]/_SumOfudsWeights;
      bookScatter2D(24, 5, 2, true)->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multPhi[2]/_SumOfcWeights;
      bookScatter2D(24, 5, 3, true)->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multPhi[3]/_SumOfbWeights;
      bookScatter2D(24, 5, 4, true)->point(0).setY(avgNumPartsBottom);
      // charm-light
      bookScatter2D(25, 5, 1, true)->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      bookScatter2D(25, 5, 2, true)->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // p
      // all
      avgNumPartsAll = _multProton[0]/sumOfWeights();
      bookScatter2D(24, 6, 1, true)->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = _multProton[1]/_SumOfudsWeights;
      bookScatter2D(24, 6, 2, true)->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multProton[2]/_SumOfcWeights;
      bookScatter2D(24, 6, 3, true)->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multProton[3]/_SumOfbWeights;
      bookScatter2D(24, 6, 4, true)->point(0).setY(avgNumPartsBottom);
      // charm-light
      bookScatter2D(25, 6, 1, true)->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      bookScatter2D(25, 6, 2, true)->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // Lambda
      // all
      avgNumPartsAll = _multLambda[0]/sumOfWeights();
      bookScatter2D(24, 7, 1, true)->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = _multLambda[1]/_SumOfudsWeights;
      bookScatter2D(24, 7, 2, true)->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = _multLambda[2]/_SumOfcWeights;
      bookScatter2D(24, 7, 3, true)->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = _multLambda[3]/_SumOfbWeights;
      bookScatter2D(24, 7, 4, true)->point(0).setY(avgNumPartsBottom);
      // charm-light
      bookScatter2D(25, 7, 1, true)->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      bookScatter2D(25, 7, 2, true)->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
    }

    //@}


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles. Used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.
    double _SumOfudsWeights, _SumOfcWeights, _SumOfbWeights;
    vector<double> _multPiPlus, _multKPlus, _multK0, _multKStar0, _multPhi, _multProton, _multLambda;

    Histo1DPtr _histXpPiPlusSig, _histXpPiPlusN;
    Histo1DPtr _histXpKPlusSig, _histXpKPlusN;
    Histo1DPtr _histXpProtonSig, _histXpProtonN;
    Histo1DPtr _histXpChargedN;
    Histo1DPtr _histXpK0N, _histXpLambdaN;
    Histo1DPtr _histXpKStar0N, _histXpPhiN;
    Histo1DPtr _histXpPiPlusLight, _histXpPiPlusCharm, _histXpPiPlusBottom;
    Histo1DPtr _histXpKPlusLight, _histXpKPlusCharm, _histXpKPlusBottom;
    Histo1DPtr _histXpKStar0Light, _histXpKStar0Charm, _histXpKStar0Bottom;
    Histo1DPtr _histXpProtonLight, _histXpProtonCharm, _histXpProtonBottom;
    Histo1DPtr _histXpLambdaLight, _histXpLambdaCharm, _histXpLambdaBottom;
    Histo1DPtr _histXpK0Light, _histXpK0Charm, _histXpK0Bottom;
    Histo1DPtr _histXpPhiLight, _histXpPhiCharm, _histXpPhiBottom;
    Histo1DPtr _tempXpKPlusCharm , _tempXpKPlusLight;
    Histo1DPtr _tempXpKStar0Charm, _tempXpKStar0Light;
    Histo1DPtr _tempXpProtonCharm, _tempXpProtonLight;
    Histo1DPtr _histRPiPlus, _histRPiMinus;
    Histo1DPtr _histRKS0, _histRKSBar0;
    Histo1DPtr _histRKPlus, _histRKMinus;
    Histo1DPtr _histRProton, _histRPBar;
    Histo1DPtr _histRLambda, _histRLBar;

    Scatter2DPtr _h_Xp_PiPl_Ch, _h_Xp_KPl_Ch,  _h_Xp_Pr_Ch;
    Scatter2DPtr _h_Xp_PiPlCh_PiPlLi, _h_Xp_PiPlBo_PiPlLi;
    Scatter2DPtr _h_Xp_KPlCh_KPlLi, _h_Xp_KPlBo_KPlLi;
    Scatter2DPtr _h_Xp_KS0Ch_KS0Li, _h_Xp_KS0Bo_KS0Li;
    Scatter2DPtr _h_Xp_PrCh_PrLi, _h_Xp_PrBo_PrLi;
    Scatter2DPtr _h_Xp_LaCh_LaLi, _h_Xp_LaBo_LaLi;
    Scatter2DPtr _h_Xp_K0Ch_K0Li, _h_Xp_K0Bo_K0Li;
    Scatter2DPtr _h_Xp_PhiCh_PhiLi, _h_Xp_PhiBo_PhiLi;

    Scatter2DPtr _h_PiM_PiP, _h_KSBar0_KS0, _h_KM_KP, _h_Pr_PBar, _h_Lam_LBar;

    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SLD_1999_S3743934);

}
