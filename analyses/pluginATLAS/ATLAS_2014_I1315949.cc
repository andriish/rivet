// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  class ATLAS_2014_I1315949 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1315949);

    void init() {

      FinalState fs;

      ZFinder zfinder(fs, Cuts::abseta<2.4 && Cuts::pT>20.0*GeV, PID::MUON, 66*GeV, 116*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder, "ZFinder");

      ChargedFinalState cfs( zfinder.remainingFinalState() );
      declare(cfs, "cfs");


      book(_h_pTsum_tow    , 67, 1, 1);
      book(_h_pTsum_trv    , 68, 1, 1);
      book(_h_pTsum_away   , 69, 1, 1);
      book(_h_pTsum_tmin   , 70, 1, 1);
      book(_h_pTsum_tmax   , 71, 1, 1);
      book(_h_pTsum_tdif   ,125, 1, 1);

      book(_h_Nchg_tow     , 72, 1, 1);
      book(_h_Nchg_trv     , 73, 1, 1);
      book(_h_Nchg_away    , 74, 1, 1);
      book(_h_Nchg_tmin    , 75, 1, 1);
      book(_h_Nchg_tmax    , 82, 1, 1);
      book(_h_Nchg_tdif    ,126, 1, 1);

      book(_h_pTavg_tow    ,113, 1, 1);
      book(_h_pTavg_trv    ,114, 1, 1);
      book(_h_pTavg_away   ,115, 1, 1);

      book(_h_pTavgvsmult_tow , 116, 1, 1);
      book(_h_pTavgvsmult_trv , 117, 1, 1);
      book(_h_pTavgvsmult_away, 118, 1, 1);


      // Book sumpt and nch histos
      for (size_t id = 0; id < 6.; ++id) {
        book(_h_ptSum_1D[0][id], 76 + id, 1, 1);
        book(_h_ptSum_1D[1][id],107 + id, 1, 1);
        book(_h_ptSum_1D[2][id],119 + id, 1, 1);
        book(_h_ptSum_1D[3][id],127 + id, 1, 1);
        book(_h_Nchg_1D[0][id],  83 + id, 1, 1);
        book(_h_Nchg_1D[1][id],  89 + id, 1, 1);
        book(_h_Nchg_1D[2][id],  95 + id, 1, 1);
        book(_h_Nchg_1D[3][id], 101 + id, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");

      if (zfinder.bosons().size() != 1) vetoEvent;

      double  Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
      double  Zphi  = zfinder.bosons()[0].momentum().phi(MINUSPI_PLUSPI);
      double  Zmass = zfinder.bosons()[0].momentum().mass()/GeV;
      if(Zmass < 66. || Zmass > 116.) vetoEvent;

      // Initialise counters for Nch and sumPt for all regions
      int nTowards(0), nTransverse(0), nLeft(0), nRight(0), nTrmin(0), nTrmax(0), nAway(0);
      double ptSumTowards(0.0), ptSumTransverse(0.0), ptSumLeft(0.0), ptSumRight(0.0),
             ptSumTrmin(0.0), ptSumTrmax(0.0), ptSumAway(0.0);

      // The charged particles
      Particles particles = apply<ChargedFinalState>(event, "cfs").particlesByPt(
          Cuts::pT > 0.5*GeV && Cuts::abseta <2.5);

      // Loop over charged particles with pT>500 MeV and |eta|<2.5
      for (const Particle& p : particles) {
        double dphi = p.momentum().phi(MINUSPI_PLUSPI) - Zphi;
        double pT   = p.momentum().pT();

        // Get multiples of 2pi right
        for(; std::fabs(dphi) > M_PI; dphi += (dphi > 0. ? -2.*M_PI : 2.*M_PI) );

        // Towards region
        if( std::fabs(dphi) < M_PI/3. ) {
          nTowards++;
          ptSumTowards += pT;
        }
        // Transverse region
        else if( std::fabs(dphi) < 2.*M_PI/3. ) {
          nTransverse++;
          ptSumTransverse += pT;
          if(dphi > 0.) {
            nRight++;
            ptSumRight += pT;
          }
          else {
            nLeft++;
            ptSumLeft += pT;
          }
        }
        // Away region
        else {
          nAway++;
          ptSumAway += pT;
        }
      }

      // TransMAX, TransMIN regions
      if (ptSumLeft > ptSumRight) {
        ptSumTrmax = ptSumLeft;
        ptSumTrmin = ptSumRight;
        nTrmax     = nLeft;
        nTrmin     = nRight;
      }
      else {
	ptSumTrmax = ptSumRight;
	ptSumTrmin = ptSumLeft;
	nTrmax     = nRight;
	nTrmin     = nLeft;
      }

      // min max regions have difference are than all other regions
      const double area = 5.*2./3.*M_PI;

      // Fill sumPt vs. Zpt region profiles
      _h_pTsum_tow->fill( Zpt, ptSumTowards/area);
      _h_pTsum_trv->fill( Zpt, ptSumTransverse/area);
      _h_pTsum_away->fill(Zpt, ptSumAway/area);
      _h_pTsum_tmin->fill(Zpt, ptSumTrmin/(0.5*area));
      _h_pTsum_tmax->fill(Zpt, ptSumTrmax/(0.5*area));
      _h_pTsum_tdif->fill(Zpt, (ptSumTrmax - ptSumTrmin)/(0.5*area));

      // Fill Nch vs. Zpt region profiles
      _h_Nchg_tow->fill( Zpt, nTowards/area);
      _h_Nchg_trv->fill( Zpt, nTransverse/area);
      _h_Nchg_away->fill(Zpt, nAway/area);
      _h_Nchg_tmin->fill(Zpt, nTrmin/(0.5*area));
      _h_Nchg_tmax->fill(Zpt, nTrmax/(0.5*area));
      _h_Nchg_tdif->fill(Zpt, (nTrmax - nTrmin)/(0.5*area));


      // Fill <pT> vs. ZpT profiles
      _h_pTavg_tow->fill( Zpt, nTowards    > 0.? ptSumTowards/nTowards       : 0.);
      _h_pTavg_trv->fill( Zpt, nTransverse > 0.? ptSumTransverse/nTransverse : 0.);
      _h_pTavg_away->fill(Zpt, nAway       > 0.? ptSumAway/nAway             : 0.);

      // Fill <Nch> vs. ZpT profiles
      _h_pTavgvsmult_tow->fill( nTowards,    nTowards    > 0.? ptSumTowards/nTowards       : 0.);
      _h_pTavgvsmult_trv->fill( nTransverse, nTransverse > 0.? ptSumTransverse/nTransverse : 0.);
      _h_pTavgvsmult_away->fill(nAway,       nAway       > 0.? ptSumAway/nAway             : 0.);

      // Determine Zpt region histo to fill
      int i_bin(0);
      if (inRange(Zpt,0,5)        ) i_bin=0;
      if (inRange(Zpt,5,10)       ) i_bin=1;
      if (inRange(Zpt,10,20)      ) i_bin=2;
      if (inRange(Zpt,20,50)      ) i_bin=3;
      if (inRange(Zpt,50,110)     ) i_bin=4;
      if (Zpt>110) i_bin=5;

      // SumPt histos for Zpt region
      _h_ptSum_1D[0][i_bin]->fill(ptSumTowards/area);
      _h_ptSum_1D[1][i_bin]->fill(ptSumTransverse/area);
      _h_ptSum_1D[2][i_bin]->fill(ptSumTrmin/(0.5*area));
      _h_ptSum_1D[3][i_bin]->fill(ptSumTrmax/(0.5*area));

      // Nch histos for Zpt region
      _h_Nchg_1D[0][i_bin]->fill(nTowards/area);
      _h_Nchg_1D[1][i_bin]->fill(nTransverse/area);
      _h_Nchg_1D[2][i_bin]->fill(nTrmin/(0.5*area));
      _h_Nchg_1D[3][i_bin]->fill(nTrmax/(0.5*area));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(int i_reg = 0; i_reg < 4; i_reg++) {
        for(int i_bin = 0; i_bin < 6; i_bin++) {
          normalize( _h_ptSum_1D[i_reg][i_bin] );
          normalize( _h_Nchg_1D[ i_reg][i_bin] );
        }
      }
    }


  private:

    Profile1DPtr _h_pTsum_tow,
                 _h_pTsum_trv,
                 _h_pTsum_away,
                 _h_pTsum_tmin,
                 _h_pTsum_tmax,
                 _h_pTsum_tdif,

                 _h_Nchg_tow,
                 _h_Nchg_trv,
                 _h_Nchg_away,
                 _h_Nchg_tmin,
                 _h_Nchg_tmax,
                 _h_Nchg_tdif,

                 _h_pTavg_tow,
                 _h_pTavg_trv,
                 _h_pTavg_away,
                 _h_pTavgvsmult_tow,
                 _h_pTavgvsmult_trv,
                 _h_pTavgvsmult_away;

    Histo1DPtr   _h_ptSum_1D[4][6], _h_Nchg_1D[4][6];


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1315949);

}
