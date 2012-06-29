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
  class ARGUS_1993_S2653028 : public Analysis {
  public:

    ARGUS_1993_S2653028()
      : Analysis("ARGUS_1993_S2653028"), _weightSum(0.)
    { }


    void analyze(const Event& e) {
      const double weight = e.weight();

      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");

      // find the upsilons
      ParticleVector upsilons;
      // first in unstable final state
      foreach (const Particle& p, ufs.particles())
        if(p.pdgId()==300553) upsilons.push_back(p);
      // then in whole event if fails
      if(upsilons.empty()) {
        foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
          if(p->pdg_id()!=300553) continue;
          const GenVertex* pv = p->production_vertex();
          bool passed = true;
          if (pv) {
            for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
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

      // find an upsilons
      foreach (const Particle& p, upsilons) {
        _weightSum += weight;
        vector<GenParticle *> pionsA,pionsB,protonsA,protonsB,kaons;
        // find the decay products we want
        findDecayProducts(p.genParticle(),pionsA,pionsB,
                          protonsA,protonsB,kaons);
        LorentzTransform cms_boost;
        if(p.momentum().vector3().mod()>0.001)
          cms_boost = LorentzTransform(-p.momentum().boostVector());
        for(unsigned int ix=0;ix<pionsA.size();++ix) {
          FourMomentum ptemp(pionsA[ix]->momentum());
          FourMomentum p2 = cms_boost.transform(ptemp);
          double pcm =
            cms_boost.transform(ptemp).vector3().mod();
          _histPiA->fill(pcm,weight);
        }
        _multPiA->fill(10.58,double(pionsA.size())*weight);
        for(unsigned int ix=0;ix<pionsB.size();++ix) {
          double pcm =
            cms_boost.transform(FourMomentum(pionsB[ix]->momentum())).vector3().mod();
          _histPiB->fill(pcm,weight);
        }
        _multPiB->fill(10.58,double(pionsB.size())*weight);
        for(unsigned int ix=0;ix<protonsA.size();++ix) {
          double pcm =
            cms_boost.transform(FourMomentum(protonsA[ix]->momentum())).vector3().mod();
          _histpA->fill(pcm,weight);
        }
        _multpA->fill(10.58,double(protonsA.size())*weight);
        for(unsigned int ix=0;ix<protonsB.size();++ix) {
          double pcm =
            cms_boost.transform(FourMomentum(protonsB[ix]->momentum())).vector3().mod();
          _histpB->fill(pcm,weight);
        }
        _multpB->fill(10.58,double(protonsB.size())*weight);
        for(unsigned int ix=0;ix<kaons.size();++ix) {
          double pcm =
            cms_boost.transform(FourMomentum(kaons[ix]->momentum())).vector3().mod();
          _histKA->fill(pcm,weight);
          _histKB->fill(pcm,weight);
        }
        _multK->fill(10.58,double(kaons.size())*weight);
      }
    } // analyze

    void finalize() {
      if(_weightSum>0.) {
        scale(_histPiA, 1./_weightSum);
        scale(_histPiB, 1./_weightSum);
        scale(_histKA , 1./_weightSum);
        scale(_histKB , 1./_weightSum);
        scale(_histpA , 1./_weightSum);
        scale(_histpB , 1./_weightSum);
        scale(_multPiA, 1./_weightSum);
        scale(_multPiB, 1./_weightSum);
        scale(_multK  , 1./_weightSum);
        scale(_multpA , 1./_weightSum);
        scale(_multpB , 1./_weightSum);
      }
    } // finalize


    void init() {
      addProjection(UnstableFinalState(), "UFS");

      // spectra
      _histPiA = bookHisto1D(1, 1, 1);
      _histPiB = bookHisto1D(2, 1, 1);
      _histKA  = bookHisto1D(3, 1, 1);
      _histKB  = bookHisto1D(6, 1, 1);
      _histpA  = bookHisto1D(4, 1, 1);
      _histpB  = bookHisto1D(5, 1, 1);
      // multiplicities
      _multPiA = bookHisto1D( 7, 1, 1);
      _multPiB = bookHisto1D( 8, 1, 1);
      _multK   = bookHisto1D( 9, 1, 1);
      _multpA  = bookHisto1D(10, 1, 1);
      _multpB  = bookHisto1D(11, 1, 1);
    } // init

  private:

    //@{
    // count of weights
    double _weightSum;
    // Histograms
    // spectra
    Histo1DPtr _histPiA;
    Histo1DPtr _histPiB;
    Histo1DPtr _histKA;
    Histo1DPtr _histKB;
    Histo1DPtr _histpA;
    Histo1DPtr _histpB;
    // multiplicities
    Histo1DPtr _multPiA;
    Histo1DPtr _multPiB;
    Histo1DPtr _multK;
    Histo1DPtr _multpA;
    Histo1DPtr _multpB;
    //@}

    void findDecayProducts(const GenParticle & p,
                           vector<GenParticle *> & pionsA,
                           vector<GenParticle *> & pionsB,
                           vector<GenParticle *> & protonsA,
                           vector<GenParticle *> & protonsB,
                           vector<GenParticle *> & kaons) {
      int parentId = p.pdg_id();
      const GenVertex* dv = p.end_vertex();
      for (GenVertex::particles_out_const_iterator
             pp = dv->particles_out_const_begin();
           pp != dv->particles_out_const_end(); ++pp) {
        int id = abs((*pp)->pdg_id());
        if(id==PIPLUS) {
          if(parentId != LAMBDA && parentId != K0S) {
            pionsA.push_back(*pp);
            pionsB.push_back(*pp);
          }
          else
            pionsB.push_back(*pp);
        }
        else if(id==PROTON) {
          if(parentId != LAMBDA && parentId != K0S) {
            protonsA.push_back(*pp);
            protonsB.push_back(*pp);
          }
          else
            protonsB.push_back(*pp);
        }
        else if(id==KPLUS) {
          kaons.push_back(*pp);
        }
        else if((*pp)->end_vertex())
          findDecayProducts(**pp,pionsA,pionsB,
                            protonsA,protonsB,kaons);
      }
    }
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1993_S2653028);

}
