// -*- C++ -*-
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Babar chamronium spectra
  /// @author Peter Richardson
  class BABAR_2003_I593379 : public Analysis {
  public:

    BABAR_2003_I593379() 
      : Analysis("BABAR_2003_I593379"), _weightSum(0.)
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
	// find the charmonium resonances
	vector<GenParticle *> allJpsi,primaryJpsi,Psiprime,
	  all_chi_c1,all_chi_c2,primary_chi_c1,primary_chi_c2;
	findDecayProducts(p.genParticle(),allJpsi,primaryJpsi,Psiprime,
			  all_chi_c1,all_chi_c2,primary_chi_c1,primary_chi_c2);
	LorentzTransform cms_boost(-p.momentum().boostVector());
	for(unsigned int ix=0;ix<allJpsi.size();++ix) {
	  double pcm = 
	    cms_boost.transform(FourMomentum(allJpsi[ix]->momentum())).vector3().mod();
	  _hist_all_Jpsi->fill(pcm,weight);
	}
	_mult_JPsi->fill(10.58,weight*double(allJpsi.size()));
	for(unsigned int ix=0;ix<primaryJpsi.size();++ix) {
	  double pcm = 
	    cms_boost.transform(FourMomentum(primaryJpsi[ix]->momentum())).vector3().mod();
	  _hist_primary_Jpsi->fill(pcm,weight);
	}
	_mult_JPsi_direct->fill(10.58,weight*double(primaryJpsi.size()));
	for(unsigned int ix=0;ix<Psiprime.size();++ix) {
	  double pcm = 
	    cms_boost.transform(FourMomentum(Psiprime[ix]->momentum())).vector3().mod();
	  _hist_Psi_prime->fill(pcm,weight);
	}
	_mult_Psi2S->fill(10.58,weight*double(Psiprime.size()));
	for(unsigned int ix=0;ix<all_chi_c1.size();++ix) {
	  double pcm = 
	    cms_boost.transform(FourMomentum(all_chi_c1[ix]->momentum())).vector3().mod();
	  _hist_chi_c1->fill(pcm,weight);
	}
	_mult_chi_c1->fill(10.58,weight*double(all_chi_c1.size()));
	_mult_chi_c1_direct->fill(10.58,weight*double(primary_chi_c1.size()));
	for(unsigned int ix=0;ix<all_chi_c2.size();++ix) {
	  double pcm = 
	    cms_boost.transform(FourMomentum(all_chi_c2[ix]->momentum())).vector3().mod();
	  _hist_chi_c2->fill(pcm,weight);
	}
	_mult_chi_c2->fill(10.58,weight*double(all_chi_c2.size()));
	_mult_chi_c2_direct->fill(10.58,weight*double(primary_chi_c2.size()));
      }
    } // analyze

    void finalize() {

      scale(_hist_all_Jpsi    ,0.5*0.1/_weightSum);
      scale(_hist_chi_c1      ,0.5*0.1/_weightSum);
      scale(_hist_chi_c2      ,0.5*0.1/_weightSum);
      scale(_hist_Psi_prime   ,0.5*0.1/_weightSum);
      scale(_hist_primary_Jpsi,0.5*0.1/_weightSum);
      scale(_mult_JPsi         ,0.5*100./_weightSum);
      scale(_mult_JPsi_direct  ,0.5*100./_weightSum);
      scale(_mult_chi_c1       ,0.5*100./_weightSum);
      scale(_mult_chi_c1_direct,0.5*100./_weightSum);
      scale(_mult_chi_c2       ,0.5*100./_weightSum);
      scale(_mult_chi_c2_direct,0.5*100./_weightSum);
      scale(_mult_Psi2S        ,0.5*100./_weightSum);
    } // finalize


    void init() {
      addProjection(UnstableFinalState(), "UFS");

      _mult_JPsi          = bookHistogram1D(1, 1, 1);
      _mult_JPsi_direct   = bookHistogram1D(1, 1, 2);
      _mult_chi_c1        = bookHistogram1D(1, 1, 3);
      _mult_chi_c1_direct = bookHistogram1D(1, 1, 4);
      _mult_chi_c2        = bookHistogram1D(1, 1, 5);
      _mult_chi_c2_direct = bookHistogram1D(1, 1, 6);
      _mult_Psi2S         = bookHistogram1D(1, 1, 7);
      _hist_all_Jpsi      = bookHistogram1D(6, 1, 1);
      _hist_chi_c1        = bookHistogram1D(7, 1, 1);
      _hist_chi_c2        = bookHistogram1D(7, 1, 2);
      _hist_Psi_prime     = bookHistogram1D(8, 1, 1);
      _hist_primary_Jpsi  = bookHistogram1D(10, 1, 1);
    } // init

  private:

    //@{
    // count of weights
    double _weightSum;
    /// Histograms
    AIDA::IHistogram1D* _hist_all_Jpsi;
    AIDA::IHistogram1D* _hist_chi_c1;
    AIDA::IHistogram1D* _hist_chi_c2;
    AIDA::IHistogram1D* _hist_Psi_prime;
    AIDA::IHistogram1D* _hist_primary_Jpsi;

    AIDA::IHistogram1D* _mult_JPsi;
    AIDA::IHistogram1D* _mult_JPsi_direct;
    AIDA::IHistogram1D* _mult_chi_c1;
    AIDA::IHistogram1D* _mult_chi_c1_direct;
    AIDA::IHistogram1D* _mult_chi_c2;
    AIDA::IHistogram1D* _mult_chi_c2_direct;
    AIDA::IHistogram1D* _mult_Psi2S;
    //@}

    void findDecayProducts(const GenParticle & p,
			   vector<GenParticle *> & allJpsi,
			   vector<GenParticle *> & primaryJpsi,
			   vector<GenParticle *> & Psiprime,
			   vector<GenParticle *> & all_chi_c1,
			   vector<GenParticle *> & all_chi_c2,
			   vector<GenParticle *> & primary_chi_c1,
			   vector<GenParticle *> & primary_chi_c2) {
      const GenVertex* dv = p.end_vertex();
      bool isOnium(false);
      for (GenVertex::particles_in_const_iterator pp = dv->particles_in_const_begin() ;
	   pp != dv->particles_in_const_end() ; ++pp) {
	int id = (*pp)->pdg_id();
	id = id%1000;
	id -= id%10;
	id /= 10;
	if(id==44) isOnium = true;
      }
      for (GenVertex::particles_out_const_iterator
	     pp = dv->particles_out_const_begin();
	   pp != dv->particles_out_const_end(); ++pp) {
	int id = (*pp)->pdg_id();
	if(id==100443) {
	  Psiprime.push_back(*pp);
	}
	else if(id==20443) {
	  all_chi_c1.push_back(*pp);
	  if(!isOnium) primary_chi_c1.push_back(*pp);
	}
	else if(id==445) {
	  all_chi_c2.push_back(*pp);
	  if(!isOnium) primary_chi_c2.push_back(*pp);
	}
	else if(id==443) {
	  allJpsi.push_back(*pp);
	  if(!isOnium) primaryJpsi.push_back(*pp);
	}
	if((*pp)->end_vertex()) {
	  findDecayProducts(**pp,allJpsi,primaryJpsi,Psiprime,
			    all_chi_c1,all_chi_c2,primary_chi_c1,primary_chi_c2);
	}
      }
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2003_I593379);

}
