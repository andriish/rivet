// -*- C++ -*-
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/ParticleName.hh"

namespace Rivet {


  /// @brief BABAR tau lepton to three charged hadrons
  /// @author Peter Richardson
  class BABAR_2007_S7266081 : public Analysis {
  public:

    BABAR_2007_S7266081() 
      : Analysis("BABAR_2007_S7266081"), _weight_total(0.),
	_weight_pipippi(0.),_weight_Kpipi(0.),_weight_KpiK(0.),_weight_KKK(0.)
    { }


    void analyze(const Event& e) {

      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");

      // find the taus
      ParticleVector taus;
      foreach (const Particle& p, ufs.particles()) {
	if(abs(p.pdgId())!=15) continue;
	_weight_total += 1.;
	ParticleVector pip,pim,Kp,Km;
	unsigned int nstable = 0.;
	// get the boost to the rest frame
	LorentzTransform cms_boost;
	if(p.momentum().vector3().mod()>0.001)
	  cms_boost = LorentzTransform(-p.momentum().boostVector());
	// find the decay products we want
	findDecayProducts(p.genParticle(),nstable,pip,pim,Kp,Km);
	if(p.pdgId()<0) {
	  swap(pip,pim);
	  swap(Kp ,Km );
	}
	if(nstable!=4) continue;
	// pipipi
	if(pim.size()==2&&pip.size()==1) {
	  _weight_pipippi += 1.;
	  _hist_pipipi_pipipi->
	    fill((pip[0].momentum()+pim[0].momentum()+pim[1].momentum()).mass(),1.);
	  _hist_pipipi_pipi->
	    fill((pip[0].momentum()+pim[0].momentum()).mass(),1.);
	  _hist_pipipi_pipi->
	    fill((pip[0].momentum()+pim[1].momentum()).mass(),1.);
	}
	else if(pim.size()==1&&pip.size()==1&&Km.size()==1) {
	  _weight_Kpipi += 1.;
	  _hist_Kpipi_Kpipi->
	    fill((pim[0].momentum()+pip[0].momentum()+Km[0].momentum()).mass(),1.);
	  _hist_Kpipi_Kpi->
	    fill((pip[0].momentum()+Km[0].momentum()).mass(),1.);
	  _hist_Kpipi_pipi->
	    fill((pim[0].momentum()+pip[0].momentum()).mass(),1.);
	}
	else if(Kp.size()==1&&Km.size()==1&&pim.size()==1) {
	  _weight_KpiK += 1.;
	  _hist_KpiK_KpiK->
	    fill((Kp[0].momentum()+Km[0].momentum()+pim[0].momentum()).mass(),1.);
	  _hist_KpiK_KK->  
	    fill((Kp[0].momentum()+Km[0].momentum()).mass(),1.);
	  _hist_KpiK_piK->
	    fill((Kp[0].momentum()+pim[0].momentum()).mass(),1.);
	}
	else if(Kp.size()==1&&Km.size()==2) {
	  _weight_KKK += 1.;
	  _hist_KKK_KKK->
	    fill((Kp[0].momentum()+Km[0].momentum()+Km[1].momentum()).mass(),1.);
	  _hist_KKK_KK->
	    fill((Kp[0].momentum()+Km[0].momentum()).mass(),1.);
	  _hist_KKK_KK->
	    fill((Kp[0].momentum()+Km[1].momentum()).mass(),1.);
	}
      }
    } // analyze

    void finalize() {
      if(_weight_pipippi>0.) {
	scale(_hist_pipipi_pipipi, 1./_weight_pipippi);
	scale(_hist_pipipi_pipi  ,0.5/_weight_pipippi);
      }
      if(_weight_Kpipi>0.) {
	scale(_hist_Kpipi_Kpipi  , 1./_weight_Kpipi);
	scale(_hist_Kpipi_Kpi    , 1./_weight_Kpipi);
	scale(_hist_Kpipi_pipi   , 1./_weight_Kpipi);
      }
      if(_weight_KpiK>0.) {
	scale(_hist_KpiK_KpiK    , 1./_weight_KpiK);
	scale(_hist_KpiK_KK      , 1./_weight_KpiK);
	scale(_hist_KpiK_piK     , 1./_weight_KpiK);
      }
      if(_weight_KKK>0.) {
	scale(_hist_KKK_KKK      , 1./_weight_KKK);
	scale(_hist_KKK_KK       ,0.5/_weight_KKK);
      }
      AIDA::IDataPointSet * br_pipipi = bookDataPointSet(11,1,1);
      br_pipipi->point(0)->coordinate(1)->setValue     ( 100.*_weight_pipippi/_weight_total);
      br_pipipi->point(0)->coordinate(1)->setErrorPlus ( 100.*sqrt(_weight_pipippi)/_weight_total);
      br_pipipi->point(0)->coordinate(1)->setErrorMinus( 100.*sqrt(_weight_pipippi)/_weight_total);
      AIDA::IDataPointSet * br_Kpipi  = bookDataPointSet(12,1,1);
      br_Kpipi->point(0)->coordinate(1)->setValue     ( 100.*_weight_Kpipi/_weight_total);
      br_Kpipi->point(0)->coordinate(1)->setErrorPlus ( 100.*sqrt(_weight_Kpipi)/_weight_total);
      br_Kpipi->point(0)->coordinate(1)->setErrorMinus( 100.*sqrt(_weight_Kpipi)/_weight_total);
      AIDA::IDataPointSet * br_KpiK   = bookDataPointSet(13,1,1);
      br_KpiK->point(0)->coordinate(1)->setValue     ( 100.*_weight_KpiK/_weight_total);
      br_KpiK->point(0)->coordinate(1)->setErrorPlus ( 100.*sqrt(_weight_KpiK)/_weight_total);
      br_KpiK->point(0)->coordinate(1)->setErrorMinus( 100.*sqrt(_weight_KpiK)/_weight_total);
      AIDA::IDataPointSet * br_KKK    = bookDataPointSet(14,1,1);
      br_KKK->point(0)->coordinate(1)->setValue     ( 100.*_weight_KKK/_weight_total);
      br_KKK->point(0)->coordinate(1)->setErrorPlus ( 100.*sqrt(_weight_KKK)/_weight_total);
      br_KKK->point(0)->coordinate(1)->setErrorMinus( 100.*sqrt(_weight_KKK)/_weight_total);
    } // finalize

    void init() {
      addProjection(UnstableFinalState(), "UFS");

      _hist_pipipi_pipipi = bookHistogram1D( 1,1,1);
      _hist_pipipi_pipi   = bookHistogram1D( 2,1,1);
      _hist_Kpipi_Kpipi   = bookHistogram1D( 3,1,1);
      _hist_Kpipi_Kpi     = bookHistogram1D( 4,1,1);
      _hist_Kpipi_pipi    = bookHistogram1D( 5,1,1);
      _hist_KpiK_KpiK     = bookHistogram1D( 6,1,1);
      _hist_KpiK_KK       = bookHistogram1D( 7,1,1);
      _hist_KpiK_piK      = bookHistogram1D( 8,1,1);
      _hist_KKK_KKK       = bookHistogram1D( 9,1,1);
      _hist_KKK_KK        = bookHistogram1D(10,1,1);

    } // init

  private:

    //@{
    AIDA::IHistogram1D* _hist_pipipi_pipipi;
    AIDA::IHistogram1D* _hist_pipipi_pipi  ;
    AIDA::IHistogram1D* _hist_Kpipi_Kpipi  ;
    AIDA::IHistogram1D* _hist_Kpipi_Kpi    ;
    AIDA::IHistogram1D* _hist_Kpipi_pipi   ;
    AIDA::IHistogram1D* _hist_KpiK_KpiK    ;
    AIDA::IHistogram1D* _hist_KpiK_KK      ;
    AIDA::IHistogram1D* _hist_KpiK_piK     ;
    AIDA::IHistogram1D* _hist_KKK_KKK      ;
    AIDA::IHistogram1D* _hist_KKK_KK       ;

    // count of weights
    double _weight_total,_weight_pipippi,_weight_Kpipi,_weight_KpiK,_weight_KKK;
    //@}

    void findDecayProducts(const GenParticle & p, unsigned int & nstable,
			   ParticleVector & pip, ParticleVector & pim,
			   ParticleVector &  Kp, ParticleVector & Km) {
      const GenVertex* dv = p.end_vertex();
      for (GenVertex::particles_out_const_iterator
	     pp = dv->particles_out_const_begin();
	   pp != dv->particles_out_const_end(); ++pp) {
	int id = (*pp)->pdg_id();
 	if( id == PI0 ) 
	  ++nstable;
	else if ( id == K0S) 
	  ++nstable;
	else if(id==PIPLUS) {
	  pip.push_back(Particle(**pp));
	  ++nstable;
	}
	else if(id==PIMINUS) {
	  pim.push_back(Particle(**pp));
	  ++nstable;
	}
	else if(id==KPLUS) {
	  Kp.push_back(Particle(**pp));
	  ++nstable;
	}
	else if(id==KMINUS) {
	  Km.push_back(Particle(**pp));
	  ++nstable;
	}
	else if((*pp)->end_vertex()) {
	  findDecayProducts(**pp,nstable,pip,pim,Kp,Km);
	}
	else
	  ++nstable;
      }
    }
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2007_S7266081);

}
