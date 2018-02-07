// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief ARGUS vector meson production
  /// @author Peter Richardson
  class ARGUS_1993_S2789213 : public Analysis {
  public:

    ARGUS_1993_S2789213()
      : Analysis("ARGUS_1993_S2789213")
    { }


    void init() {
      declare(UnstableFinalState(), "UFS");

      book(_mult_cont_Omega     , 1, 1, 1);
      book(_mult_cont_Rho0      , 1, 1, 2);
      book(_mult_cont_KStar0    , 1, 1, 3);
      book(_mult_cont_KStarPlus , 1, 1, 4);
      book(_mult_cont_Phi       , 1, 1, 5);

      book(_mult_Ups1_Omega     , 2, 1, 1);
      book(_mult_Ups1_Rho0      , 2, 1, 2);
      book(_mult_Ups1_KStar0    , 2, 1, 3);
      book(_mult_Ups1_KStarPlus , 2, 1, 4);
      book(_mult_Ups1_Phi       , 2, 1, 5);

      book(_mult_Ups4_Omega     , 3, 1, 1);
      book(_mult_Ups4_Rho0      , 3, 1, 2);
      book(_mult_Ups4_KStar0    , 3, 1, 3);
      book(_mult_Ups4_KStarPlus , 3, 1, 4);
      book(_mult_Ups4_Phi       , 3, 1, 5);

      book(_hist_cont_KStarPlus , 4, 1, 1);
      book(_hist_Ups1_KStarPlus , 5, 1, 1);
      book(_hist_Ups4_KStarPlus , 6, 1, 1);

      book(_hist_cont_KStar0    , 7, 1, 1);
      book(_hist_Ups1_KStar0    , 8, 1, 1);
      book(_hist_Ups4_KStar0    , 9, 1, 1);

      book(_hist_cont_Rho0      ,10, 1, 1);
      book(_hist_Ups1_Rho0      ,11, 1, 1);
      book(_hist_Ups4_Rho0      ,12, 1, 1);

      book(_hist_cont_Omega     ,13, 1, 1);
      book(_hist_Ups1_Omega     ,14, 1, 1);


      book(_weightSum_cont,"TMP/weightSumcont");
      book(_weightSum_Ups1,"TMP/weightSumUps1");
      book(_weightSum_Ups4,"TMP/weightSumUps4");
    }


    void analyze(const Event& e) {
      // Find the upsilons
      Particles upsilons;
      // First in unstable final state
      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");
      foreach (const Particle& p, ufs.particles())
        if (p.pid() == 300553 || p.pid() == 553) upsilons.push_back(p);
      // Then in whole event if that failed
      if (upsilons.empty()) {
        foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
          if (p->pdg_id() != 300553 && p->pdg_id() != 553) continue;
          const GenVertex* pv = p->production_vertex();
          bool passed = true;
          if (pv) {
            foreach (const GenParticle* pp, particles_in(pv)) {
              if ( p->pdg_id() == pp->pdg_id() ) {
                passed = false;
                break;
              }
            }
          }
          if (passed) upsilons.push_back(Particle(*p));
        }
      }

      if (upsilons.empty()) { // continuum

        _weightSum_cont->fill();
        unsigned int nOmega(0), nRho0(0), nKStar0(0), nKStarPlus(0), nPhi(0);
        foreach (const Particle& p, ufs.particles()) {
          int id = p.abspid();
          double xp = 2.*p.E()/sqrtS();
          double beta = p.p3().mod()/p.E();
          if (id == 113) {
            _hist_cont_Rho0->fill(xp, 1.0/beta);
            ++nRho0;
          }
          else if (id == 313) {
            _hist_cont_KStar0->fill(xp, 1.0/beta);
            ++nKStar0;
          }
          else if (id == 223) {
            _hist_cont_Omega->fill(xp, 1.0/beta);
            ++nOmega;
          }
          else if (id == 323) {
            _hist_cont_KStarPlus->fill(xp,1.0/beta);
            ++nKStarPlus;
          }
          else if (id == 333) {
            ++nPhi;
          }
        }
        /// @todo Replace with Counters and fill one-point Scatters at the end
        _mult_cont_Omega    ->fill(10.45, nOmega    );
        _mult_cont_Rho0     ->fill(10.45, nRho0     );
        _mult_cont_KStar0   ->fill(10.45, nKStar0   );
        _mult_cont_KStarPlus->fill(10.45, nKStarPlus);
        _mult_cont_Phi      ->fill(10.45, nPhi      );

      } else { // found an upsilon

        foreach (const Particle& ups, upsilons) {
          const int parentId = ups.pid();
          (parentId == 553 ? _weightSum_Ups1 : _weightSum_Ups4)->fill();
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups.genParticle(),unstable);
          /// @todo Update to new LT mk* functions
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 0.001)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          double mass = ups.mass();
          unsigned int nOmega(0),nRho0(0),nKStar0(0),nKStarPlus(0),nPhi(0);
          foreach(const Particle & p , unstable) {
            int id = p.abspid();
            FourMomentum p2 = cms_boost.transform(p.momentum());
            double xp = 2.*p2.E()/mass;
            double beta = p2.p3().mod()/p2.E();
            if (id == 113) {
              if (parentId == 553) _hist_Ups1_Rho0->fill(xp,1.0/beta);
              else                 _hist_Ups4_Rho0->fill(xp,1.0/beta);
              ++nRho0;
            }
            else if (id == 313) {
              if (parentId == 553) _hist_Ups1_KStar0->fill(xp,1.0/beta);
              else                 _hist_Ups4_KStar0->fill(xp,1.0/beta);
              ++nKStar0;
            }
            else if (id == 223) {
              if (parentId == 553) _hist_Ups1_Omega->fill(xp,1.0/beta);
              ++nOmega;
            }
            else if (id == 323) {
              if (parentId == 553) _hist_Ups1_KStarPlus->fill(xp,1.0/beta);
              else                 _hist_Ups4_KStarPlus->fill(xp,1.0/beta);
              ++nKStarPlus;
            }
            else if (id == 333) {
              ++nPhi;
            }
          }
          if (parentId == 553) {
            _mult_Ups1_Omega    ->fill(9.46,nOmega    );
            _mult_Ups1_Rho0     ->fill(9.46,nRho0     );
            _mult_Ups1_KStar0   ->fill(9.46,nKStar0   );
            _mult_Ups1_KStarPlus->fill(9.46,nKStarPlus);
            _mult_Ups1_Phi      ->fill(9.46,nPhi      );
          }
          else {
            _mult_Ups4_Omega    ->fill(10.58,nOmega    );
            _mult_Ups4_Rho0     ->fill(10.58,nRho0     );
            _mult_Ups4_KStar0   ->fill(10.58,nKStar0   );
            _mult_Ups4_KStarPlus->fill(10.58,nKStarPlus);
            _mult_Ups4_Phi      ->fill(10.58,nPhi      );
          }
        }
      }

    }


    void finalize() {
      if (_weightSum_cont->val() > 0.) {
        /// @todo Replace with Counters and fill one-point Scatters at the end
        scale(_mult_cont_Omega    , 1. / *_weightSum_cont);
        scale(_mult_cont_Rho0     , 1. / *_weightSum_cont);
        scale(_mult_cont_KStar0   , 1. / *_weightSum_cont);
        scale(_mult_cont_KStarPlus, 1. / *_weightSum_cont);
        scale(_mult_cont_Phi      , 1. / *_weightSum_cont);
        scale(_hist_cont_KStarPlus, 1. / *_weightSum_cont);
        scale(_hist_cont_KStar0   , 1. / *_weightSum_cont);
        scale(_hist_cont_Rho0     , 1. / *_weightSum_cont);
        scale(_hist_cont_Omega    , 1. / *_weightSum_cont);
      }
      if (_weightSum_Ups1->val() > 0.) {
        /// @todo Replace with Counters and fill one-point Scatters at the end
        scale(_mult_Ups1_Omega    , 1. / *_weightSum_Ups1);
        scale(_mult_Ups1_Rho0     , 1. / *_weightSum_Ups1);
        scale(_mult_Ups1_KStar0   , 1. / *_weightSum_Ups1);
        scale(_mult_Ups1_KStarPlus, 1. / *_weightSum_Ups1);
        scale(_mult_Ups1_Phi      , 1. / *_weightSum_Ups1);
        scale(_hist_Ups1_KStarPlus, 1. / *_weightSum_Ups1);
        scale(_hist_Ups1_KStar0   , 1. / *_weightSum_Ups1);
        scale(_hist_Ups1_Rho0     , 1. / *_weightSum_Ups1);
        scale(_hist_Ups1_Omega    , 1. / *_weightSum_Ups1);
      }
      if (_weightSum_Ups4->val() > 0.) {
        /// @todo Replace with Counters and fill one-point Scatters at the end
        scale(_mult_Ups4_Omega    , 1. / *_weightSum_Ups4);
        scale(_mult_Ups4_Rho0     , 1. / *_weightSum_Ups4);
        scale(_mult_Ups4_KStar0   , 1. / *_weightSum_Ups4);
        scale(_mult_Ups4_KStarPlus, 1. / *_weightSum_Ups4);
        scale(_mult_Ups4_Phi      , 1. / *_weightSum_Ups4);
        scale(_hist_Ups4_KStarPlus, 1. / *_weightSum_Ups4);
        scale(_hist_Ups4_KStar0   , 1. / *_weightSum_Ups4);
        scale(_hist_Ups4_Rho0     , 1. / *_weightSum_Ups4);
      }
    }


  private:

    //@{
    Histo1DPtr _mult_cont_Omega, _mult_cont_Rho0, _mult_cont_KStar0, _mult_cont_KStarPlus, _mult_cont_Phi;
    Histo1DPtr _mult_Ups1_Omega, _mult_Ups1_Rho0, _mult_Ups1_KStar0, _mult_Ups1_KStarPlus, _mult_Ups1_Phi;
    Histo1DPtr _mult_Ups4_Omega, _mult_Ups4_Rho0, _mult_Ups4_KStar0, _mult_Ups4_KStarPlus, _mult_Ups4_Phi;
    Histo1DPtr _hist_cont_KStarPlus, _hist_Ups1_KStarPlus, _hist_Ups4_KStarPlus;
    Histo1DPtr _hist_cont_KStar0, _hist_Ups1_KStar0, _hist_Ups4_KStar0   ;
    Histo1DPtr _hist_cont_Rho0, _hist_Ups1_Rho0,  _hist_Ups4_Rho0;
    Histo1DPtr _hist_cont_Omega, _hist_Ups1_Omega;

    CounterPtr _weightSum_cont,_weightSum_Ups1,_weightSum_Ups4;
    //@}


    void findDecayProducts(const GenParticle* p, Particles& unstable) {
      const GenVertex* dv = p->end_vertex();
      /// @todo Use better looping
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        int id = abs((*pp)->pdg_id());
        if (id == 113 || id == 313 || id == 323 ||
            id == 333 || id == 223 ) {
          unstable.push_back(Particle(*pp));
        }
        else if ((*pp)->end_vertex())
          findDecayProducts(*pp, unstable);
      }
    }


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1993_S2789213);

}
