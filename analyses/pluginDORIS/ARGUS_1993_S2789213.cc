// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief ARGUS vector meson production
  /// @author Peter Richardson
  class ARGUS_1993_S2789213 : public Analysis {
  public:

    ARGUS_1993_S2789213()
      : Analysis("ARGUS_1993_S2789213")
    { }


    void init() {
      declare(UnstableParticles(), "UFS");
<<<<<<< local

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
=======
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<5;++iy) {
	  ostringstream title;
	  title << "/TMP/MULT_" << ix << "_" << iy;
	  _mult[ix][iy] = bookCounter(title.str());
	}
      }
>>>>>>> graft

      book(_hist_cont_KStar0    , 7, 1, 1);
      book(_hist_Ups1_KStar0    , 8, 1, 1);
      book(_hist_Ups4_KStar0    , 9, 1, 1);

      book(_hist_cont_Rho0      ,10, 1, 1);
      book(_hist_Ups1_Rho0      ,11, 1, 1);
      book(_hist_Ups4_Rho0      ,12, 1, 1);

      book(_hist_cont_Omega     ,13, 1, 1);
      book(_hist_Ups1_Omega     ,14, 1, 1);


<<<<<<< local
      book(_weightSum_cont,"TMP/weightSumcont");
      book(_weightSum_Ups1,"TMP/weightSumUps1");
      book(_weightSum_Ups4,"TMP/weightSumUps4");
=======
      _hist_cont_Omega     = bookHisto1D(13, 1, 1);
      _hist_Ups1_Omega     = bookHisto1D(14, 1, 1);

>>>>>>> graft
    }


    void analyze(const Event& e) {
      // Find the upsilons
      // First in unstable final state
<<<<<<< local
      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");
      for (const Particle& p : ufs.particles())
        if (p.pid() == 300553 || p.pid() == 553) upsilons.push_back(p);
      // Then in whole event if that failed
      if (upsilons.empty()) {
        for(ConstGenParticlePtr p: HepMCUtils::particles(e.genEvent())) {
          if (p->pdg_id() != 300553 && p->pdg_id() != 553) continue;
          ConstGenVertexPtr pv = p->production_vertex();
          bool passed = true;
          if (pv) {
            for(ConstGenParticlePtr pp: HepMCUtils::particles(pv, Relatives::PARENTS)){
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
        for (const Particle& p : ufs.particles()) {
=======
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 || Cuts::pid==300553);
      // continuum
      if (upsilons.empty()) { 
        _weightSum_cont += weight;
        foreach (const Particle& p, ufs.particles()) {
>>>>>>> graft
          int id = p.abspid();
          double xp = 2.*p.E()/sqrtS();
          double beta = p.p3().mod()/p.E();
          if (id == 113) {
<<<<<<< local
            _hist_cont_Rho0->fill(xp, 1.0/beta);
            ++nRho0;
=======
            _hist_cont_Rho0->fill(xp, weight/beta);
	    _mult[0][1]->fill(weight);
>>>>>>> graft
          }
          else if (id == 313) {
<<<<<<< local
            _hist_cont_KStar0->fill(xp, 1.0/beta);
            ++nKStar0;
=======
            _hist_cont_KStar0->fill(xp, weight/beta);
	    _mult[0][2]->fill(weight);
>>>>>>> graft
          }
          else if (id == 223) {
<<<<<<< local
            _hist_cont_Omega->fill(xp, 1.0/beta);
            ++nOmega;
=======
            _hist_cont_Omega->fill(xp, weight/beta);
	    _mult[0][0]->fill(weight);
>>>>>>> graft
          }
          else if (id == 323) {
<<<<<<< local
            _hist_cont_KStarPlus->fill(xp,1.0/beta);
            ++nKStarPlus;
=======
            _hist_cont_KStarPlus->fill(xp,weight/beta);
	    _mult[0][3]->fill(weight);
>>>>>>> graft
          }
          else if (id == 333) {
	    _mult[0][4]->fill(weight);
          }
        }
<<<<<<< local
        /// @todo Replace with Counters and fill one-point Scatters at the end
        _mult_cont_Omega    ->fill(10.45, nOmega    );
        _mult_cont_Rho0     ->fill(10.45, nRho0     );
        _mult_cont_KStar0   ->fill(10.45, nKStar0   );
        _mult_cont_KStarPlus->fill(10.45, nKStarPlus);
        _mult_cont_Phi      ->fill(10.45, nPhi      );

      } else { // found an upsilon

        for (const Particle& ups : upsilons) {
=======
      }
      // found an upsilon
      else {
        foreach (const Particle& ups, upsilons) {
>>>>>>> graft
          const int parentId = ups.pid();
          (parentId == 553 ? _weightSum_Ups1 : _weightSum_Ups4)->fill();
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups,unstable);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 0.001)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          double mass = ups.mass();
<<<<<<< local
          unsigned int nOmega(0),nRho0(0),nKStar0(0),nKStarPlus(0),nPhi(0);
          for(const Particle & p : unstable) {
=======
          foreach(const Particle & p , unstable) {
>>>>>>> graft
            int id = p.abspid();
            FourMomentum p2 = cms_boost.transform(p.momentum());
            double xp = 2.*p2.E()/mass;
            double beta = p2.p3().mod()/p2.E();
            if (id == 113) {
<<<<<<< local
              if (parentId == 553) _hist_Ups1_Rho0->fill(xp,1.0/beta);
              else                 _hist_Ups4_Rho0->fill(xp,1.0/beta);
              ++nRho0;
=======
              if (parentId == 553) {
		_hist_Ups1_Rho0->fill(xp,weight/beta);
		_mult[1][1]->fill(weight);
	      }
              else {
		_hist_Ups4_Rho0->fill(xp,weight/beta);
		_mult[2][1]->fill(weight);
	      }
>>>>>>> graft
            }
            else if (id == 313) {
<<<<<<< local
              if (parentId == 553) _hist_Ups1_KStar0->fill(xp,1.0/beta);
              else                 _hist_Ups4_KStar0->fill(xp,1.0/beta);
              ++nKStar0;
=======
              if (parentId == 553) {
		_hist_Ups1_KStar0->fill(xp,weight/beta);
		_mult[1][2]->fill(weight);
	      }
              else {
		_hist_Ups4_KStar0->fill(xp,weight/beta);
		_mult[2][2]->fill(weight);
	      }
>>>>>>> graft
            }
            else if (id == 223) {
<<<<<<< local
              if (parentId == 553) _hist_Ups1_Omega->fill(xp,1.0/beta);
              ++nOmega;
=======
              if (parentId == 553) {
		_hist_Ups1_Omega->fill(xp,weight/beta);
		_mult[1][0]->fill(weight);
	      }
              else {
		_mult[2][0]->fill(weight);
	      }
>>>>>>> graft
            }
            else if (id == 323) {
<<<<<<< local
              if (parentId == 553) _hist_Ups1_KStarPlus->fill(xp,1.0/beta);
              else                 _hist_Ups4_KStarPlus->fill(xp,1.0/beta);
              ++nKStarPlus;
=======
              if (parentId == 553) {
		_hist_Ups1_KStarPlus->fill(xp,weight/beta);
		_mult[1][3]->fill(weight);
	      }
              else {
		_hist_Ups4_KStarPlus->fill(xp,weight/beta);
		_mult[2][3]->fill(weight);
	      }
>>>>>>> graft
            }
            else if (id == 333) {
              if (parentId == 553) {
		_mult[1][4]->fill(weight);
	      }
	      else {
		_mult[2][4]->fill(weight);
	      }
            }
          }
<<<<<<< local
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
=======
>>>>>>> graft
        }
      }

    }


    void finalize() {
<<<<<<< local
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
=======
      // multiplicities
      vector<double> scales = {_weightSum_cont,_weightSum_Ups1,_weightSum_Ups4};
      for(unsigned int ix=0;ix<3;++ix) {
	if(scales[ix]<=0.) continue;
	for(unsigned int iy=0;iy<5;++iy) {
	  Scatter2DPtr scatter = bookScatter2D(ix+1, 1, iy+1, true);
	  scale(_mult[ix][iy],1./scales[ix]);
	  scatter->point(0).setY(_mult[ix][iy]->val(),_mult[ix][iy]->err());
	}
      }
      // spectra
      if (_weightSum_cont > 0.) {
        scale(_hist_cont_KStarPlus, 1./_weightSum_cont);
        scale(_hist_cont_KStar0   , 1./_weightSum_cont);
        scale(_hist_cont_Rho0     , 1./_weightSum_cont);
        scale(_hist_cont_Omega    , 1./_weightSum_cont);
>>>>>>> graft
      }
<<<<<<< local
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
=======
      if (_weightSum_Ups1 > 0.) {
        scale(_hist_Ups1_KStarPlus, 1./_weightSum_Ups1);
        scale(_hist_Ups1_KStar0   , 1./_weightSum_Ups1);
        scale(_hist_Ups1_Rho0     , 1./_weightSum_Ups1);
        scale(_hist_Ups1_Omega    , 1./_weightSum_Ups1);
>>>>>>> graft
      }
<<<<<<< local
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
=======
      if (_weightSum_Ups4 > 0.) {
        scale(_hist_Ups4_KStarPlus, 1./_weightSum_Ups4);
        scale(_hist_Ups4_KStar0   , 1./_weightSum_Ups4);
        scale(_hist_Ups4_Rho0     , 1./_weightSum_Ups4);
>>>>>>> graft
      }
    }


  private:

    //@{
    Histo1DPtr _hist_cont_KStarPlus, _hist_Ups1_KStarPlus, _hist_Ups4_KStarPlus;
    Histo1DPtr _hist_cont_KStar0, _hist_Ups1_KStar0, _hist_Ups4_KStar0   ;
    Histo1DPtr _hist_cont_Rho0, _hist_Ups1_Rho0,  _hist_Ups4_Rho0;
    Histo1DPtr _hist_cont_Omega, _hist_Ups1_Omega;

<<<<<<< local
    CounterPtr _weightSum_cont,_weightSum_Ups1,_weightSum_Ups4;
=======
    CounterPtr _mult[3][5];



    
    double _weightSum_cont,_weightSum_Ups1,_weightSum_Ups4;
>>>>>>> graft
    //@}


   /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pdgId());
	if (id == 113 || id == 313 || id == 323 ||
            id == 333 || id == 223 ) {
	  unstable.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1993_S2789213);

}
