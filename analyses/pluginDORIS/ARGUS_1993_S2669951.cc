// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Production of the $\eta'(958)$ and $f_0(980)$ in $e^+e^-$ annihilation in the Upsilon region
  /// @author Peter Richardson
  class ARGUS_1993_S2669951 : public Analysis {
  public:

    ARGUS_1993_S2669951()
      : Analysis("ARGUS_1993_S2669951"),
        _weightSum_cont(0.),
        _weightSum_Ups1(0.),
        _weightSum_Ups2(0.)
    {   }


    void init() {
      declare(UnstableParticles(), "UFS");

      // histograms
      _hist_cont_f0 = bookHisto1D(2, 1, 1);
      _hist_Ups1_f0 = bookHisto1D(3, 1, 1);
      _hist_Ups2_f0 = bookHisto1D(4, 1, 1);

      // counters
      _count_etaPrime_highZ.push_back(bookCounter("/TMP/EtaHighCont"));
      _count_etaPrime_highZ.push_back(bookCounter("/TMP/EtaHighUps1"));
      _count_etaPrime_allZ.push_back(bookCounter("/TMP/EtaAllCont"));
      _count_etaPrime_allZ.push_back(bookCounter("/TMP/EtaAllUps1"));
      _count_etaPrime_allZ.push_back(bookCounter("/TMP/EtaAllUps2"));
      _count_f0.push_back(bookCounter("/TMP/f0Cont"));
      _count_f0.push_back(bookCounter("/TMP/f0Ups1"));
      _count_f0.push_back(bookCounter("/TMP/f0Ups2"));
    }


    void analyze(const Event& e) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 or Cuts::pid==100553);
      // Finding done, now fill counters
      const double weight = e.weight();
      // Continuum
      if (upsilons.empty()) { 
        MSG_DEBUG("No Upsilons found => continuum event");

        _weightSum_cont += weight;
        foreach (const Particle& p, ufs.particles()) {
          const int id = p.pdgId();
          const double xp = 2.*p.E()/sqrtS();
          const double beta = p.p3().mod() / p.E();
          if (id == 9010221) {
            _hist_cont_f0->fill(xp, weight/beta);
	    _count_f0[2]->fill(weight);
          } else if (id == 331) {
            if (xp > 0.35) _count_etaPrime_highZ[1]->fill(weight);
	    _count_etaPrime_allZ[2]->fill(weight);
          }
        }
      }
      // Upsilon(s) found
      else { 
        MSG_DEBUG("Upsilons found => resonance event");

        for (const Particle& ups : upsilons) {
          const int parentId = ups.pid();
	  if(parentId==553) {
	    _weightSum_Ups1 += weight;
	  }
	  else {
	    _weightSum_Ups2 += weight;
	  }
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups, unstable);
	  // boost to rest frame (if required)
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();
	  // loop over decay products
          foreach(const Particle& p, unstable) {
            const int id = p.pdgId();
            const FourMomentum p2 = cms_boost.transform(p.momentum());
            const double xp = 2.*p2.E()/mass;
            const double beta = p2.p3().mod()/p2.E();
            if (id == 9010221) {
	      if(parentId == 553 ) {
		_hist_Ups1_f0->fill(xp, weight/beta);
		_count_f0[0]->fill(weight);
	      }
	      else {
		_hist_Ups2_f0->fill(xp, weight/beta);
		_count_f0[1]->fill(weight);
	      }
	    }
	    else if ( id == 331 ) {
	      if (parentId == 553) {
		if (xp > 0.35) _count_etaPrime_highZ[0]->fill(weight);
		_count_etaPrime_allZ[0]->fill(weight);
	      }
	      else {
		_count_etaPrime_allZ[1]->fill(weight);
	      }
	    }
	  }
	}
      }
    }


    void finalize() {
      // High-Z eta' multiplicity
      Scatter2DPtr s111 = bookScatter2D(1, 1, 1, true);
      // Point at 9.460
      if (_weightSum_Ups1 > 0) {
	scale(_count_etaPrime_highZ[0] , 1./_weightSum_Ups1);
        s111->point(0).setY(_count_etaPrime_highZ[0]->val(),
			    _count_etaPrime_highZ[0]->err());
      }
      // Point at 9.905
      if (_weightSum_cont > 0) {
	scale(_count_etaPrime_highZ[1] , 1./_weightSum_cont);
	s111->point(1).setY(_count_etaPrime_highZ[1]->val(),
			    _count_etaPrime_highZ[1]->err());
      }

      // All-Z eta' multiplicity
      Scatter2DPtr s112 = bookScatter2D(1, 1, 2, true);
      // Point at 9.460
      if (_weightSum_Ups1 > 0) {
	scale(_count_etaPrime_allZ[0] , 1./ _weightSum_Ups1);
        s112->point(0).setY(_count_etaPrime_allZ[0]->val(),
			    _count_etaPrime_allZ[0]->err());
      }
      // Point at 9.905
      if (_weightSum_cont > 0) {
	scale(_count_etaPrime_allZ[2] , 1. / _weightSum_cont);
	s112->point(1).setY(_count_etaPrime_allZ[2]->val(),
			    _count_etaPrime_allZ[2]->err());
      }
      // Point at 10.02
      if (_weightSum_Ups2 > 0) {
	scale(_count_etaPrime_allZ[1] , 1. / _weightSum_Ups2);
	s112->point(2).setY(_count_etaPrime_allZ[1]->val(),
			    _count_etaPrime_allZ[1]->err());
      }
      
      // f0 multiplicity
      Scatter2DPtr s511 = bookScatter2D(5, 1, 1, true);
      // Point at 9.46
      if (_weightSum_Ups1 > 0) {
	scale(_count_f0[0] , 1./ _weightSum_Ups1);
        s511->point(0).setY(_count_f0[0]->val(),
			    _count_f0[0]->err());
      }
      // Point at 10.02
      if (_weightSum_Ups2 > 0) {
	scale(_count_f0[1] , 1./ _weightSum_Ups2);
        s511->point(1).setY(_count_f0[1]->val(),
			    _count_f0[1]->err());
      }
      // Point at 10.45
      if (_weightSum_cont > 0) {
	scale(_count_f0[2] , 1./ _weightSum_cont);
        s511->point(2).setY(_count_f0[2]->val(),
			    _count_f0[2]->err());
      }
      // Scale histos
      if (_weightSum_cont > 0.) scale(_hist_cont_f0, 1./_weightSum_cont);
      if (_weightSum_Ups1 > 0.) scale(_hist_Ups1_f0, 1./_weightSum_Ups1);
      if (_weightSum_Ups2 > 0.) scale(_hist_Ups2_f0, 1./_weightSum_Ups2);
    }


  private:

    /// @name Counters
    //@{
    vector<CounterPtr> _count_etaPrime_highZ, _count_etaPrime_allZ, _count_f0;
    double _weightSum_cont,_weightSum_Ups1,_weightSum_Ups2;
    //@}


    /// Histos
    Histo1DPtr _hist_cont_f0, _hist_Ups1_f0, _hist_Ups2_f0;


    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = p.pdgId();
	if (id == 331 || id == 9010221) {
	  unstable.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1993_S2669951);

}
