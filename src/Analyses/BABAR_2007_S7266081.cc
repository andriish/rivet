// -*- C++ -*-
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/RivetYODA.hh"
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

      // Find the taus
      Particles taus;
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");
      foreach (const Particle& p, ufs.particles()) {
        if(abs(p.pdgId())!=15) continue;
        _weight_total += 1.;
        Particles pip,pim,Kp,Km;
        unsigned int nstable = 0.;
        // get the boost to the rest frame
        LorentzTransform cms_boost;
        if(p.momentum().vector3().mod()>0.001)
          cms_boost = LorentzTransform(-p.momentum().boostVector());
        // find the decay products we want
        findDecayProducts(p.genParticle(), nstable, pip, pim, Kp, Km);
        if (p.pdgId()<0) {
          swap(pip,pim);
          swap(Kp ,Km );
        }
        if (nstable!=4) continue;
        // pipipi
        if (pim.size()==2&&pip.size()==1) {
          _weight_pipippi += 1.;
          _hist_pipipi_pipipi->
            fill((pip[0].momentum()+pim[0].momentum()+pim[1].momentum()).mass(),1.);
          _hist_pipipi_pipi->
            fill((pip[0].momentum()+pim[0].momentum()).mass(),1.);
          _hist_pipipi_pipi->
            fill((pip[0].momentum()+pim[1].momentum()).mass(),1.);
        }
        else if (pim.size()==1&&pip.size()==1&&Km.size()==1) {
          _weight_Kpipi += 1.;
          _hist_Kpipi_Kpipi->
            fill((pim[0].momentum()+pip[0].momentum()+Km[0].momentum()).mass(),1.);
          _hist_Kpipi_Kpi->
            fill((pip[0].momentum()+Km[0].momentum()).mass(),1.);
          _hist_Kpipi_pipi->
            fill((pim[0].momentum()+pip[0].momentum()).mass(),1.);
        }
        else if (Kp.size()==1&&Km.size()==1&&pim.size()==1) {
          _weight_KpiK += 1.;
          _hist_KpiK_KpiK->
            fill((Kp[0].momentum()+Km[0].momentum()+pim[0].momentum()).mass(),1.);
          _hist_KpiK_KK->
            fill((Kp[0].momentum()+Km[0].momentum()).mass(),1.);
          _hist_KpiK_piK->
            fill((Kp[0].momentum()+pim[0].momentum()).mass(),1.);
        }
        else if (Kp.size()==1&&Km.size()==2) {
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
      // @todo YODA
      //AIDA::IDataPointSet * br_pipipi = bookDataPointSet(11,1,1);
      //br_pipipi->point(0)->coordinate(1)->setValue     ( 100.*_weight_pipippi/_weight_total);
      //br_pipipi->point(0)->coordinate(1)->setErrorPlus ( 100.*sqrt(_weight_pipippi)/_weight_total);
      //br_pipipi->point(0)->coordinate(1)->setErrorMinus( 100.*sqrt(_weight_pipippi)/_weight_total);
      //AIDA::IDataPointSet * br_Kpipi  = bookDataPointSet(12,1,1);
      //br_Kpipi->point(0)->coordinate(1)->setValue     ( 100.*_weight_Kpipi/_weight_total);
      //br_Kpipi->point(0)->coordinate(1)->setErrorPlus ( 100.*sqrt(_weight_Kpipi)/_weight_total);
      //br_Kpipi->point(0)->coordinate(1)->setErrorMinus( 100.*sqrt(_weight_Kpipi)/_weight_total);
      //AIDA::IDataPointSet * br_KpiK   = bookDataPointSet(13,1,1);
      //br_KpiK->point(0)->coordinate(1)->setValue     ( 100.*_weight_KpiK/_weight_total);
      //br_KpiK->point(0)->coordinate(1)->setErrorPlus ( 100.*sqrt(_weight_KpiK)/_weight_total);
      //br_KpiK->point(0)->coordinate(1)->setErrorMinus( 100.*sqrt(_weight_KpiK)/_weight_total);
      //AIDA::IDataPointSet * br_KKK    = bookDataPointSet(14,1,1);
      //br_KKK->point(0)->coordinate(1)->setValue     ( 100.*_weight_KKK/_weight_total);
      //br_KKK->point(0)->coordinate(1)->setErrorPlus ( 100.*sqrt(_weight_KKK)/_weight_total);
      //br_KKK->point(0)->coordinate(1)->setErrorMinus( 100.*sqrt(_weight_KKK)/_weight_total);
    } // finalize


    void init() {
      addProjection(UnstableFinalState(), "UFS");

      _hist_pipipi_pipipi = bookHisto1D( 1,1,1);
      _hist_pipipi_pipi   = bookHisto1D( 2,1,1);
      _hist_Kpipi_Kpipi   = bookHisto1D( 3,1,1);
      _hist_Kpipi_Kpi     = bookHisto1D( 4,1,1);
      _hist_Kpipi_pipi    = bookHisto1D( 5,1,1);
      _hist_KpiK_KpiK     = bookHisto1D( 6,1,1);
      _hist_KpiK_KK       = bookHisto1D( 7,1,1);
      _hist_KpiK_piK      = bookHisto1D( 8,1,1);
      _hist_KKK_KKK       = bookHisto1D( 9,1,1);
      _hist_KKK_KK        = bookHisto1D(10,1,1);

    } // init


  private:

    //@{
    Histo1DPtr _hist_pipipi_pipipi;
    Histo1DPtr _hist_pipipi_pipi  ;
    Histo1DPtr _hist_Kpipi_Kpipi  ;
    Histo1DPtr _hist_Kpipi_Kpi    ;
    Histo1DPtr _hist_Kpipi_pipi   ;
    Histo1DPtr _hist_KpiK_KpiK    ;
    Histo1DPtr _hist_KpiK_KK      ;
    Histo1DPtr _hist_KpiK_piK     ;
    Histo1DPtr _hist_KKK_KKK      ;
    Histo1DPtr _hist_KKK_KK       ;

    // count of weights
    double _weight_total,_weight_pipippi,_weight_Kpipi,_weight_KpiK,_weight_KKK;
    //@}

    void findDecayProducts(const GenParticle* p,
                           unsigned int & nstable,
                           Particles& pip, Particles& pim,
                           Particles&  Kp, Particles& Km) {
      const GenVertex* dv = p->end_vertex();
      /// @todo Use better looping
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
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
          findDecayProducts(*pp, nstable, pip, pim, Kp, Km);
        }
        else
          ++nstable;
      }
    }


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2007_S7266081);

}
