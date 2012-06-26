// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "LWH/AIDataPointSet.h"
#include "LWH/Histogram1D.h"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include <iostream>
#include <fstream>

namespace Rivet {


  class ATLAS_2012_I1091481 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1091481()
      : Analysis("ATLAS_2012_I1091481")
    {
    }


  public:

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here
      ChargedFinalState cfs(-2.5, 2.5, 0.1*GeV);
      addProjection(cfs,"CFS");
      ChargedFinalState cfs5(-2.5, 2.5, 0.5*GeV);
      addProjection(cfs5,"CFS5");

      // pion mass;
      pim = 0.1396;

      // collision energy
      isqrts = -1;
      if (fuzzyEquals(sqrtS(), 900*GeV)) isqrts = 2;
      else if (fuzzyEquals(sqrtS(), 7*TeV)) isqrts = 1;
      assert(isqrts >= 0);

      // graphs
      _se_incl = bookDataPointSet(isqrts, 1, 1);
      _se_lpte = bookDataPointSet(isqrts, 1, 2);
      _se_lptd = bookDataPointSet(isqrts, 1, 3);
      _seta_incl = bookDataPointSet(isqrts, 2, 1);
      _seta_lpte = bookDataPointSet(isqrts, 2, 2);
      _seta_lptd = bookDataPointSet(isqrts, 2, 3);

      nbe = _se_incl->size();
      nbeta = _seta_incl->size();

      // parameters
      x_e.resize(nbe);
      x_eta.resize(nbeta);

      sec.resize(nbe);
      ses.resize(nbe);
      setac.resize(nbeta);
      setas.resize(nbeta);


      float bw = (binEdges(isqrts,1,1).back()-binEdges(isqrts,1,1).front())/(nbe-1);
      float lb = binEdges(isqrts,1,1).front()-0.5*bw;
      float ub = binEdges(isqrts,1,1).back()+0.5*bw;
      for (int i=0; i< nbe; i++) x_e[i]=lb+0.5*(2*i+1)*bw;

      // histograms
      h_se_incl = new LWH::Histogram1D(nbe,lb,ub);
      h_se_lpte = new LWH::Histogram1D(nbe,lb,ub);
      h_se_lptd = new LWH::Histogram1D(nbe,lb,ub);


      bw = (binEdges(isqrts,2,1).back()-binEdges(isqrts,2,1).front())/(nbeta-1);
      lb = binEdges(isqrts,2,1).front()-0.5*bw;
      ub = binEdges(isqrts,2,1).back()+0.5*bw;
      for (int j=0; j< nbeta; j++) x_eta[j]=lb+0.5*(2*j+1)*bw;

      // histograms
      h_seta_incl = new LWH::Histogram1D(nbeta,lb,ub);
      h_seta_lpte = new LWH::Histogram1D(nbeta,lb,ub);
      h_seta_lptd = new LWH::Histogram1D(nbeta,lb,ub);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& had = applyProjection<ChargedFinalState>(event, "CFS");

      vector<HepMC::GenParticle> hdrs;

      double ptmax=0.;
      foreach (const Particle& p, had.particles()) {
        const GenParticle& gp=p.genParticle();
        if (gp.momentum().perp()>ptmax) ptmax=gp.momentum().perp();
        // eta ordering
        if (!hdrs.size()) hdrs.push_back(gp);
        else {
          vector<HepMC::GenParticle>::iterator pit=hdrs.begin();
          while (pit!=hdrs.end() && (*pit).momentum().eta()>gp.momentum().eta()) pit++;
          hdrs.insert(pit,gp);
        }
      }

      //for (int i=0; i<hdrs.size(); i++) {
      //        cout <<i << ","<<hdrs[i].momentum().perp()<<","<<hdrs[i].momentum().eta() << endl;
      //}

      if (ptmax>10. || hdrs.size()<11 ) return;

      // inclusive+low pt enhanced
      for (int i=0;i<nbe;i++) {
        sec[i]=0.;
        ses[i]=0.;
      }
      for (int i=0;i<nbeta;i++) {
        setac[i]=0.;
        setas[i]=0.;
      }

      int nch = hdrs.size();

      double evis=0.;
      for (unsigned int i=0;i<hdrs.size();i++) {
        const GenParticle& gp=hdrs[i];
        double ap = gp.momentum().perp()/sin(gp.momentum().theta());
        double en = sqrt(ap*ap+pim*pim);
        evis += 0.5*en;
        double ph = gp.momentum().phi();
        double eta = gp.momentum().eta();
        for (int ii=0;ii<nbe;ii++) {
          sec[ii]+=cos(x_e[ii]*evis-ph);
          ses[ii]+=sin(x_e[ii]*evis-ph);
        }
        for (int jj=0;jj<nbeta;jj++) {
          setac[jj]+=cos(x_eta[jj]*eta-ph);
          setas[jj]+=sin(x_eta[jj]*eta-ph);
        }
        evis += 0.5*en;
      }

      for (int i=0;i<nbe;i++) {
        h_se_incl->fill(x_e[i],(sec[i]*sec[i]+ses[i]*ses[i])/nch-1.);
        if (ptmax < 1.) h_se_lpte->fill(x_e[i],(sec[i]*sec[i]+ses[i]*ses[i])/nch-1.);
      }
      for (int j=0;j<nbeta;j++) {
        h_seta_incl->fill(x_eta[j],(setac[j]*setac[j]+setas[j]*setas[j])/nch-1.);
        if (ptmax < 1.) h_seta_lpte->fill(x_eta[j],(setac[j]*setac[j]+setas[j]*setas[j])/nch-1.);
      }

      // low-pT depleted

      const ChargedFinalState& had5 = applyProjection<ChargedFinalState>(event, "CFS5");

      vector<HepMC::GenParticle> hdrs5;

      ptmax=0.;
      foreach (const Particle& p, had5.particles()) {
        const GenParticle& gp=p.genParticle();
        if (gp.momentum().perp()>ptmax) ptmax=gp.momentum().perp();
        // eta ordering
        if (!hdrs5.size()) hdrs5.push_back(gp);
        else {
          vector<HepMC::GenParticle>::iterator pit=hdrs5.begin();
          while (pit!=hdrs5.end() && (*pit).momentum().eta()>gp.momentum().eta()) pit++;
          hdrs5.insert(pit,gp);
        }
      }

      if (ptmax>10. || hdrs5.size()<11 ) return;

      for (int i=0;i<nbe;i++) {
        sec[i]=0.;
        ses[i]=0.;
      }
      for (int i=0;i<nbeta;i++) {
        setac[i]=0.;
        setas[i]=0.;
      }

      nch = hdrs5.size();

      evis=0.;
      for (unsigned int i=0;i<hdrs5.size();i++) {
        const GenParticle& gp=hdrs5[i];
        double ap = gp.momentum().perp()/sin(gp.momentum().theta());
        double en = sqrt(ap*ap+pim*pim);
        evis += 0.5*en;
        double ph = gp.momentum().phi();
        double eta = gp.momentum().eta();
        for (int ii=0;ii<nbe;ii++) {
          sec[ii]+=cos(x_e[ii]*evis-ph);
          ses[ii]+=sin(x_e[ii]*evis-ph);
        }
        for (int jj=0;jj<nbeta;jj++) {
          setac[jj]+=cos(x_eta[jj]*eta-ph);
          setas[jj]+=sin(x_eta[jj]*eta-ph);
        }
        evis += 0.5*en;
      }

      for (int i=0;i<nbe;i++) {
        h_se_lptd->fill(x_e[i],(sec[i]*sec[i]+ses[i]*ses[i])/nch-1.);
      }
      for (int j=0;j<nbeta;j++) {
        h_seta_lptd->fill(x_eta[j],(setac[j]*setac[j]+setas[j]*setas[j])/nch-1.);
      }

      return;
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // number of events :
      int nev_incl=h_se_incl->entries()/nbe;
      int nev_lpte=h_se_lpte->entries()/nbe;
      int nev_lptd=h_se_lptd->entries()/nbe;

      // inclusive
      for (int i=0; i<nbe; i++) {
        sec[i]= h_se_incl->binHeight(i)/nev_incl;
        ses[i]= h_se_incl->binError(i)/nev_incl;
      }
      _se_incl->setCoordinate(1,sec,ses);

      for (int i=0; i<nbeta; i++){
        setac[i]= h_seta_incl->binHeight(i)/nev_incl;
        setas[i]= h_seta_incl->binError(i)/nev_incl;
      }
      _seta_incl->setCoordinate(1,setac,setas);

      // low-pT enhanced
      for (int i=0; i<nbe; i++) {
        sec[i]= h_se_lpte->binHeight(i)/nev_lpte;
        ses[i]= h_se_lpte->binError(i)/nev_lpte;
      }
      _se_lpte->setCoordinate(1,sec,ses);

      for (int i=0; i<nbeta; i++){
        setac[i]= h_seta_lpte->binHeight(i)/nev_lpte;
        setas[i]= h_seta_lpte->binError(i)/nev_lpte;
      }
      _seta_lpte->setCoordinate(1,setac,setas);

      // low-pT depleted
      for (int i=0; i<nbe; i++) {
        sec[i]= h_se_lptd->binHeight(i)/nev_lptd;
        ses[i]= h_se_lptd->binError(i)/nev_lptd;
      }
      _se_lptd->setCoordinate(1,sec,ses);

      for (int i=0; i<nbeta; i++){
        setac[i]= h_seta_lptd->binHeight(i)/nev_lptd;
        setas[i]= h_seta_lptd->binError(i)/nev_lptd;
      }
      _seta_lptd->setCoordinate(1,setac,setas);

      // save event counts ?

    }

    //@}


  private:

    int isqrts;
    int nbe;
    int nbeta;

    AIDA::IDataPointSet*       _se_incl;
    AIDA::IDataPointSet*       _se_lpte;
    AIDA::IDataPointSet*       _se_lptd;
    AIDA::IDataPointSet*       _seta_incl;
    AIDA::IDataPointSet*       _seta_lpte;
    AIDA::IDataPointSet*       _seta_lptd;

    LWH::Histogram1D* h_se_incl;
    LWH::Histogram1D* h_se_lpte;
    LWH::Histogram1D* h_se_lptd;
    LWH::Histogram1D* h_seta_incl;
    LWH::Histogram1D* h_seta_lpte;
    LWH::Histogram1D* h_seta_lptd;

    double pim;

    // arrays
    std::vector<float>         x_e;
    std::vector<float>         x_eta;

    std::vector<double>        sec;
    std::vector<double>        ses;
    std::vector<double>        setac;
    std::vector<double>        setas;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1091481);


}
