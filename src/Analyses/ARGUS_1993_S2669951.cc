// -*- C++ -*-
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/RivetYODA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/ParticleName.hh"

namespace Rivet {


  /// @brief BELLE pi+/-, K+/- and proton/antiproton spectrum at Upsilon(4S)
  /// @author Peter Richardson
  class ARGUS_1993_S2669951 : public Analysis {
  public:

    ARGUS_1993_S2669951()
      : Analysis("ARGUS_1993_S2669951"), _count_etaPrime_highZ(2,0.),
        _count_etaPrime_allZ(3,0.), _count_f0(3,0.),
        _weightSum_cont(0.),_weightSum_Ups1(0.),_weightSum_Ups2(0.)
    { }


    void analyze(const Event& e) {
      const double weight = e.weight();

      const Beam beamproj = applyProjection<Beam>(e, "Beams");
      const double s = sqr(beamproj.sqrtS());
      const double roots = sqrt(s);
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");

      // find the upsilons
      Particles upsilons;
      // first in unstable final state
      foreach (const Particle& p, ufs.particles())
        if (p.pdgId()==553 || p.pdgId()==100553 ) upsilons.push_back(p);
      // then in whole event if fails
      if (upsilons.empty()) {
        foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
          if( p->pdg_id() != 553 && p->pdg_id() != 100553 ) continue;
          const GenVertex* pv = p->production_vertex();
          bool passed = true;
          if (pv) {
            for (GenVertex::particles_in_const_iterator
                   pp = pv->particles_in_const_begin() ;
                 pp != pv->particles_in_const_end() ; ++pp) {
              if ( p->pdg_id() == (*pp)->pdg_id() ) {
                passed = false;
                break;
              }
            }
          }
          if(passed) upsilons.push_back(Particle(*p));
        }
      }
      // continuum
      if(upsilons.empty()) {
        _weightSum_cont += weight;
        unsigned int nEtaA(0),nEtaB(0),nf0(0);
        foreach (const Particle& p, ufs.particles()) {
          int id = abs(p.pdgId());
          double xp = 2.*p.momentum().t()/roots;
          double beta = p.momentum().vector3().mod()/p.momentum().t();
          if(id==9010221) {
            _hist_cont_f0->fill(xp,weight/beta);
            ++nf0;
          }
          else if(id==331) {
            if(xp>0.35) ++nEtaA;
            ++nEtaB;
          }
        }
        _count_f0[2]             += nf0*weight;
        _count_etaPrime_highZ[1] += nEtaA*weight;
        _count_etaPrime_allZ[2]  += nEtaB*weight;
      }
      else {
        // find an upsilons
        foreach (const Particle& ups, upsilons) {
          int parentId = ups.pdgId();
          if(parentId==553)
            _weightSum_Ups1 += weight;
          else
            _weightSum_Ups2 += weight;
          Particles unstable;
          // find the decay products we want
          findDecayProducts(ups.genParticle(), unstable);
          LorentzTransform cms_boost;
          if(ups.momentum().vector3().mod()>0.001)
            cms_boost = LorentzTransform(-ups.momentum().boostVector());
          double mass = ups.momentum().mass();
          unsigned int nEtaA(0),nEtaB(0),nf0(0);
          foreach(const Particle & p , unstable) {
            int id = abs(p.pdgId());
            FourMomentum p2 = cms_boost.transform(p.momentum());
            double xp = 2.*p2.t()/mass;
            double beta = p2.vector3().mod()/p2.t();
            if (id==9010221) {
              if(parentId==553) _hist_Ups1_f0->fill(xp,weight/beta);
              else              _hist_Ups2_f0->fill(xp,weight/beta);
              ++nf0;
            }
            else if (id==331) {
              if(xp>0.35) ++nEtaA;
              ++nEtaB;
            }
          }
          if (parentId==553) {
            _count_f0[0]             +=   nf0*weight;
            _count_etaPrime_highZ[0] += nEtaA*weight;
            _count_etaPrime_allZ[0]  += nEtaB*weight;
          }
          else {
            _count_f0[1] += nf0*weight;
            _count_etaPrime_allZ[1] += nEtaB*weight;
          }
        }
      }
    } // analyze

    void finalize() {

      // @todo YODA
      //AIDA::IDataPointSet * mult_etaPrime_highZ = bookDataPointSet( 1,1,1);
      //for (int i = 0; i < mult_etaPrime_highZ->size(); ++i) {
      //  if ( fuzzyEquals( 9.905, mult_etaPrime_highZ->point(i)->coordinate(0)->value(), 0.001) &&
      //       _weightSum_cont>0.)
      //    mult_etaPrime_highZ->point(i)->coordinate(1)->setValue( _count_etaPrime_highZ[1]/_weightSum_cont);
      //  else if ( fuzzyEquals( 9.46, mult_etaPrime_highZ->point(i)->coordinate(0)->value(), 0.001) &&
      //            _weightSum_Ups1>0.)
      //    mult_etaPrime_highZ->point(i)->coordinate(1)->setValue(_count_etaPrime_highZ[0]/_weightSum_Ups1);
      //}
      //AIDA::IDataPointSet * mult_etaPrime_allZ = bookDataPointSet( 1,1,2);
      //for (int i = 0; i < mult_etaPrime_allZ->size(); ++i) {
      //  if ( fuzzyEquals( 9.905, mult_etaPrime_allZ->point(i)->coordinate(0)->value(), 0.001) &&
      //       _weightSum_cont>0.) {
      //    mult_etaPrime_allZ->point(i)->coordinate(1)->setValue( _count_etaPrime_allZ[2]/_weightSum_cont);
      //  }
      //  else if ( fuzzyEquals( 9.46, mult_etaPrime_allZ->point(i)->coordinate(0)->value(), 0.001) &&
      //            _weightSum_Ups1>0.) {
      //    mult_etaPrime_allZ->point(i)->coordinate(1)->setValue(_count_etaPrime_allZ[0]/_weightSum_Ups1);
      //  }
      //  else if ( fuzzyEquals( 10.02, mult_etaPrime_allZ->point(i)->coordinate(0)->value(), 0.001) &&
      //            _weightSum_Ups2>0.) {
      //    mult_etaPrime_allZ->point(i)->coordinate(1)->setValue(_count_etaPrime_allZ[1]/_weightSum_Ups2);
      //  }
      //}
      //AIDA::IDataPointSet * mult_f0 = bookDataPointSet( 5,1,1);
      //for (int i = 0; i < mult_f0->size(); ++i) {
      //  if ( fuzzyEquals( 10.45, mult_f0->point(i)->coordinate(0)->value(), 0.001) &&
      //       _weightSum_cont>0.) {
      //    mult_f0->point(i)->coordinate(1)->setValue( _count_f0[2]/_weightSum_cont);
      //  }
      //  else if ( fuzzyEquals( 9.46, mult_f0->point(i)->coordinate(0)->value(), 0.001) &&
      //            _weightSum_Ups1>0.) {
      //    mult_f0->point(i)->coordinate(1)->setValue(_count_f0[0]/_weightSum_Ups1);
      //  }
      //  else if ( fuzzyEquals( 10.02, mult_f0->point(i)->coordinate(0)->value(), 0.001) &&
      //            _weightSum_Ups2>0.) {
      //    mult_f0->point(i)->coordinate(1)->setValue(_count_f0[1]/_weightSum_Ups2);
      //  }
      //}

      if(_weightSum_cont>0.) scale(_hist_cont_f0, 1./_weightSum_cont);
      if(_weightSum_Ups1>0.) scale(_hist_Ups1_f0, 1./_weightSum_Ups1);
      if(_weightSum_Ups2>0.) scale(_hist_Ups2_f0, 1./_weightSum_Ups2);

    } // finalize


    void init() {
      addProjection(Beam(), "Beams");
      addProjection(UnstableFinalState(), "UFS");

      _hist_cont_f0 = bookHisto1D( 2,1,1);
      _hist_Ups1_f0 = bookHisto1D( 3,1,1);
      _hist_Ups2_f0 = bookHisto1D( 4,1,1);
    }


  private:

    //@{
    vector<double> _count_etaPrime_highZ;
    vector<double> _count_etaPrime_allZ;
    vector<double> _count_f0;

    Histo1DPtr _hist_cont_f0;
    Histo1DPtr _hist_Ups1_f0;
    Histo1DPtr _hist_Ups2_f0;

    double _weightSum_cont,_weightSum_Ups1,_weightSum_Ups2;
    //@}


    void findDecayProducts(const GenParticle* p,
                           Particles & unstable) {
      const GenVertex* dv = p->end_vertex();
      /// @todo Use better looping
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        int id = abs((*pp)->pdg_id());
        if (id == 331 || id == 9010221) {
          unstable.push_back(Particle(*pp));
        } else if ((*pp)->end_vertex())
          findDecayProducts(*pp, unstable);
      }
    }

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1993_S2669951);

}

