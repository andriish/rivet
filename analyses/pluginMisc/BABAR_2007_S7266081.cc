// -*- C++ -*-
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief BABAR tau lepton to three charged hadrons
  /// @author Peter Richardson
  class BABAR_2007_S7266081 : public Analysis {
  public:

    BABAR_2007_S7266081()
      : Analysis("BABAR_2007_S7266081")
    {   }


    void init() {
      declare(UnstableFinalState(), "UFS");
      book(_hist_pipipi_pipipi , 1, 1, 1);
      book(_hist_pipipi_pipi   , 2, 1, 1);
      book(_hist_Kpipi_Kpipi   , 3, 1, 1);
      book(_hist_Kpipi_Kpi     , 4, 1, 1);
      book(_hist_Kpipi_pipi    , 5, 1, 1);
      book(_hist_KpiK_KpiK     , 6, 1, 1);
      book(_hist_KpiK_KK       , 7, 1, 1);
      book(_hist_KpiK_piK      , 8, 1, 1);
      book(_hist_KKK_KKK       , 9, 1, 1);
      book(_hist_KKK_KK        ,10, 1, 1);

      book(_weight_total, "weight_total");
      book(_weight_pipipi, "weight_pipipi");
      book(_weight_Kpipi, "weight_Kpipi");
      book(_weight_KpiK, "weight_KpiK");
      book(_weight_KKK, "weight_KKK");

    }


    void analyze(const Event& e) {
      // Find the taus
      Particles taus;
      foreach(const Particle& p, apply<UnstableFinalState>(e, "UFS").particles(Cuts::pid==PID::TAU)) {
        _weight_total->fill();
        Particles pip, pim, Kp, Km;
        unsigned int nstable = 0;
        // Get the boost to the rest frame
        LorentzTransform cms_boost;
        if (p.p3().mod() > 1*MeV)
          cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        // Find the decay products we want
        findDecayProducts(p.genParticle(), nstable, pip, pim, Kp, Km);
        if (p.pid() < 0) {
          swap(pip, pim);
          swap(Kp, Km );
        }
        if (nstable != 4) continue;
        // pipipi
        if (pim.size() == 2 && pip.size() == 1) {
          _weight_pipipi->fill();
          _hist_pipipi_pipipi->
            fill((pip[0].momentum()+pim[0].momentum()+pim[1].momentum()).mass());
          _hist_pipipi_pipi->
            fill((pip[0].momentum()+pim[0].momentum()).mass());
          _hist_pipipi_pipi->
            fill((pip[0].momentum()+pim[1].momentum()).mass());
        }
        else if (pim.size() == 1 && pip.size() == 1 && Km.size() == 1) {
          _weight_Kpipi->fill();
          _hist_Kpipi_Kpipi->
            fill((pim[0].momentum()+pip[0].momentum()+Km[0].momentum()).mass());
          _hist_Kpipi_Kpi->
            fill((pip[0].momentum()+Km[0].momentum()).mass());
          _hist_Kpipi_pipi->
            fill((pim[0].momentum()+pip[0].momentum()).mass());
        }
        else if (Kp.size() == 1 && Km.size() == 1 && pim.size() == 1) {
          _weight_KpiK->fill();
          _hist_KpiK_KpiK->
            fill((Kp[0].momentum()+Km[0].momentum()+pim[0].momentum()).mass());
          _hist_KpiK_KK->
            fill((Kp[0].momentum()+Km[0].momentum()).mass());
          _hist_KpiK_piK->
            fill((Kp[0].momentum()+pim[0].momentum()).mass());
        }
        else if (Kp.size() == 1 && Km.size() == 2) {
          _weight_KKK->fill();
          _hist_KKK_KKK->
            fill((Kp[0].momentum()+Km[0].momentum()+Km[1].momentum()).mass());
          _hist_KKK_KK->
            fill((Kp[0].momentum()+Km[0].momentum()).mass());
          _hist_KKK_KK->
            fill((Kp[0].momentum()+Km[1].momentum()).mass());
        }
      }
    }


    void finalize() {
      if (_weight_pipipi > 0.) {
        scale(_hist_pipipi_pipipi, 1.0/_weight_pipipi);
        scale(_hist_pipipi_pipi  , 0.5/_weight_pipipi);
      }
      if (_weight_Kpipi > 0.) {
        scale(_hist_Kpipi_Kpipi  , 1.0/_weight_Kpipi);
        scale(_hist_Kpipi_Kpi    , 1.0/_weight_Kpipi);
        scale(_hist_Kpipi_pipi   , 1.0/_weight_Kpipi);
      }
      if (_weight_KpiK > 0.) {
        scale(_hist_KpiK_KpiK    , 1.0/_weight_KpiK);
        scale(_hist_KpiK_KK      , 1.0/_weight_KpiK);
        scale(_hist_KpiK_piK     , 1.0/_weight_KpiK);
      }
      if (_weight_KKK > 0.) {
        scale(_hist_KKK_KKK      , 1.0/_weight_KKK);
        scale(_hist_KKK_KK       , 0.5/_weight_KKK);
      }
      /// @note Using autobooking for these scatters since their x values are not really obtainable from the MC data
      Scatter2DPtr tmp11, tmp12, tmp13, tmp14;
      book(tmp11, 11, 1, 1, true);
      book(tmp12, 12, 1, 1, true);
      book(tmp13, 13, 1, 1, true);
      book(tmp14, 14, 1, 1, true);      
      tmp11->point(0).setY(100*_weight_pipipi/_weight_total, 100*sqrt(double(_weight_pipipi))/_weight_total);
      tmp12->point(0).setY(100*_weight_Kpipi/_weight_total, 100*sqrt(double(_weight_Kpipi))/_weight_total);
      tmp13->point(0).setY(100*_weight_KpiK/_weight_total, 100*sqrt(double(_weight_KpiK))/_weight_total);
      tmp14->point(0).setY(100*_weight_KKK/_weight_total, 100*sqrt(double(_weight_KKK))/_weight_total);
    }


  private:

    //@{

    // Histograms
    Histo1DPtr _hist_pipipi_pipipi, _hist_pipipi_pipi;
    Histo1DPtr _hist_Kpipi_Kpipi, _hist_Kpipi_Kpi, _hist_Kpipi_pipi;
    Histo1DPtr _hist_KpiK_KpiK, _hist_KpiK_KK, _hist_KpiK_piK;
    Histo1DPtr _hist_KKK_KKK, _hist_KKK_KK;

    // Weights counters
    CounterPtr _weight_total, _weight_pipipi, _weight_Kpipi, _weight_KpiK, _weight_KKK;

    //@}

    void findDecayProducts(const GenParticle* p,
                           unsigned int & nstable,
                           Particles& pip, Particles& pim,
                           Particles& Kp, Particles& Km) {
      const GenVertex* dv = p->end_vertex();
      /// @todo Use better looping
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        int id = (*pp)->pdg_id();
        if (id == PID::PI0 )
          ++nstable;
        else if (id == PID::K0S)
          ++nstable;
        else if (id == PID::PIPLUS) {
          pip.push_back(Particle(**pp));
          ++nstable;
        }
        else if (id == PID::PIMINUS) {
          pim.push_back(Particle(**pp));
          ++nstable;
        }
        else if (id == PID::KPLUS) {
          Kp.push_back(Particle(**pp));
          ++nstable;
        }
        else if (id == PID::KMINUS) {
          Km.push_back(Particle(**pp));
          ++nstable;
        }
        else if ((*pp)->end_vertex()) {
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
