// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Production of the $\eta'(958)$ and $f_0(980)$ in $e^+e^-$ annihilation in the Upsilon region
  /// @author Peter Richardson
  class ARGUS_1993_S2669951 : public Analysis {
  public:

    ARGUS_1993_S2669951()
      : Analysis("ARGUS_1993_S2669951")
    {   }


    void init() {
      declare(UnstableFinalState(), "UFS");

      book(_weightSum_cont, "weightSum_cont");
      book(_weightSum_Ups1, "weightSum_Ups1");
      book(_weightSum_Ups2, "weightSum_Ups2");

      for ( auto i : {0,1,2} ) {
        if ( i < 2 )
          book(_count_etaPrime_highZ[i], "count_etaPrime_highz_" + to_str(i));
        book(_count_etaPrime_allZ[i], "count_etaPrime_allz_" + to_str(i));
        book(_count_f0[i], "count_f0_" + to_str(i));
      }

      book(_hist_cont_f0 ,2, 1, 1);
      book(_hist_Ups1_f0 ,3, 1, 1);
      book(_hist_Ups2_f0 ,4, 1, 1);

      book(s111, 1, 1, 1, true);
      book(s112, 1, 1, 2, true);
      book(s511, 5, 1, 1, true);

    }


    void analyze(const Event& e) {

      // Find the Upsilons among the unstables
      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");
      Particles upsilons;

      // First in unstable final state
      foreach (const Particle& p, ufs.particles())
        if (p.pid() == 553 || p.pid() == 100553)
          upsilons.push_back(p);
      // Then in whole event if fails
      if (upsilons.empty()) {
        /// @todo Replace HepMC digging with Particle::descendents etc. calls
        foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
          if ( p->pdg_id() != 553 && p->pdg_id() != 100553 ) continue;
          // Discard it if its parent has the same PDG ID code (avoid duplicates)
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


      // Finding done, now fill counters
      if (upsilons.empty()) { // Continuum
        MSG_DEBUG("No Upsilons found => continuum event");

        _weightSum_cont->fill();
        unsigned int nEtaA(0), nEtaB(0), nf0(0);
        foreach (const Particle& p, ufs.particles()) {
          const int id = p.abspid();
          const double xp = 2.*p.E()/sqrtS();
          const double beta = p.p3().mod() / p.E();
          if (id == 9010221) {
            _hist_cont_f0->fill(xp, 1.0/beta);
            nf0 += 1;
          } else if (id == 331) {
            if (xp > 0.35) nEtaA += 1;
            nEtaB += 1;
          }
        }
        _count_f0[2]            ->fill(nf0);
        _count_etaPrime_highZ[1]->fill(nEtaA);
        _count_etaPrime_allZ[2] ->fill(nEtaB);

      } else { // Upsilon(s) found
        MSG_DEBUG("Upsilons found => resonance event");

        foreach (const Particle& ups, upsilons) {
          const int parentId = ups.pid();
          ((parentId == 553) ? _weightSum_Ups1 : _weightSum_Ups2)->fill();
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups.genParticle(), unstable);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();
          unsigned int nEtaA(0), nEtaB(0), nf0(0);
          foreach(const Particle& p, unstable) {
            const int id = p.abspid();
            const FourMomentum p2 = cms_boost.transform(p.momentum());
            const double xp = 2.*p2.E()/mass;
            const double beta = p2.p3().mod()/p2.E();
            if (id == 9010221) { //< ?
              ((parentId == 553) ? _hist_Ups1_f0 : _hist_Ups2_f0)->fill(xp, 1.0/beta);
              nf0 += 1;
            } else if (id == 331) { //< ?
              if (xp > 0.35) nEtaA += 1;
              nEtaB += 1;
            }
          }
          if (parentId == 553) {
            _count_f0[0]            ->fill(  nf0);
            _count_etaPrime_highZ[0]->fill(nEtaA);
            _count_etaPrime_allZ[0] ->fill(nEtaB);
          } else {
            _count_f0[1]->fill(nf0);
            _count_etaPrime_allZ[1] ->fill(nEtaB);
          }
        }
      }
    }


    void finalize() {

      // High-Z eta' multiplicity
      if (_weightSum_Ups1 > 0) // Point at 9.460
        s111->point(0).setY(_count_etaPrime_highZ[0] / _weightSum_Ups1, 0);
      if (_weightSum_cont > 0) // Point at 9.905
        s111->point(1).setY(_count_etaPrime_highZ[1] / _weightSum_cont, 0);

      // All-Z eta' multiplicity
      if (_weightSum_Ups1 > 0) // Point at 9.460
        s112->point(0).setY(_count_etaPrime_allZ[0] / _weightSum_Ups1, 0);
      if (_weightSum_cont > 0) // Point at 9.905
        s112->point(1).setY(_count_etaPrime_allZ[2] / _weightSum_cont, 0);
      if (_weightSum_Ups2 > 0) // Point at 10.02
        s112->point(2).setY(_count_etaPrime_allZ[1] / _weightSum_Ups2, 0);

      // f0 multiplicity
      if (_weightSum_Ups1 > 0) // Point at 9.46
        s511->point(0).setY(_count_f0[0] / _weightSum_Ups1, 0);
      if (_weightSum_Ups2 > 0) // Point at 10.02
        s511->point(1).setY(_count_f0[1] / _weightSum_Ups2, 0);
      if (_weightSum_cont > 0) // Point at 10.45
        s511->point(2).setY(_count_f0[2] / _weightSum_cont, 0);

      // Scale histos
      if (_weightSum_cont > 0.) scale(_hist_cont_f0, 1./_weightSum_cont);
      if (_weightSum_Ups1 > 0.) scale(_hist_Ups1_f0, 1./_weightSum_Ups1);
      if (_weightSum_Ups2 > 0.) scale(_hist_Ups2_f0, 1./_weightSum_Ups2);
    }


  private:

    /// @name Counters
    //@{
    array<CounterPtr,3> _count_etaPrime_highZ, _count_etaPrime_allZ, _count_f0;
    CounterPtr _weightSum_cont,_weightSum_Ups1,_weightSum_Ups2;
    //@}

    Scatter2DPtr s111, s112, s511;

    /// Histos
    Histo1DPtr _hist_cont_f0, _hist_Ups1_f0, _hist_Ups2_f0;


    /// Recursively walk the HepMC tree to find decay products of @a p
    void findDecayProducts(const GenParticle* p, Particles& unstable) {
      const GenVertex* dv = p->end_vertex();
      /// @todo Use better looping
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        const int id = abs((*pp)->pdg_id());
        if (id == 331 || id == 9010221) unstable.push_back(Particle(*pp));
        else if ((*pp)->end_vertex()) findDecayProducts(*pp, unstable);
      }
    }


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1993_S2669951);

}
