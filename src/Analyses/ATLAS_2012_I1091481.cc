// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "LWH/AIDataPointSet.h"
#include "LWH/Histogram1D.h"
#include "LWH/Profile1D.h"
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
      ChargedFinalState cfs100(-2.5, 2.5, 0.1*GeV);
      addProjection(cfs100,"CFS100");
      ChargedFinalState cfs500(-2.5, 2.5, 0.5*GeV);
      addProjection(cfs500,"CFS500");


      // collision energy
      int isqrts = -1;
      if (fuzzyEquals(sqrtS(), 900*GeV)) isqrts = 2;       
      else if (fuzzyEquals(sqrtS(), 7*TeV)) isqrts = 1;       
      assert(isqrts >= 0);

      _sE_10_100   = bookHistogram1D(isqrts, 1, 1);
      _sE_1_100    = bookHistogram1D(isqrts, 1, 2);
      _sE_10_500   = bookHistogram1D(isqrts, 1, 3);

      _sEta_10_100 = bookHistogram1D(isqrts, 2, 1);
      _sEta_1_100  = bookHistogram1D(isqrts, 2, 2);
      _sEta_10_500 = bookHistogram1D(isqrts, 2, 3);
    }
   
    std::vector<double> getXj(const ParticleVector& part) {
      // Iterate over particles to get vector X_j (energy thingy with PION mass for all particles)
      //
      // X_j = 0.5*E_j + sum_{k=0}^{k<j}(E_k)
      //
      // pion mass;
      double m_pi = 0.1396;     
      
      std::vector<double> xj;
      std::vector<double> Ej;
      foreach (const Particle& p, part) {
        double pT = p.momentum().pT();
        double theta = p.momentum().theta();
        double E_j = sqrt(pow(m_pi,2) + pow(pT/sin(theta), 2));
        
        Ej.push_back(E_j);
        double temp = 0.5*E_j;
        if (xj.size()==0) xj.push_back(temp);
        else {
          for (unsigned int k=0; k< xj.size(); k++) {
            temp += Ej[k];
          }
          xj.push_back(temp);
        }
      }
      return xj;
    }
   
    // Return the stuff that needs to get filled, dependent on those parameters xi and omega
    double getSeta(const ParticleVector& part, double xi) {
      // The cores of the double sum
      double c_eta = 0.0;
     
      // This is the inner double sum
      for (unsigned int i=0; i<part.size(); i++) {
        //for (unsigned int j=0; j<part.size(); j++) {
        for (unsigned int j=0; j<i; j++) {
          if (i!=j) {
            const Particle& p_i = part[i];
            const Particle& p_j = part[j];
            double dphi = deltaPhi(p_i, p_j);
            double deta = p_i.momentum().eta() - p_j.momentum().eta();
            c_eta += cos(xi*deta    - dphi);
          }
        }
      }
      return c_eta/part.size();
    }
    
    double getSE(const ParticleVector& part, std::vector<double> Xj, double omega) {
      // The cores of the double sum
      double c_E   = 0.0;
     
      // This is the inner double sum
      for (unsigned int i=0; i<part.size(); i++) {
        //for (unsigned int j=0; j<part.size(); j++) {
        for (unsigned int j=0; j<i; j++) {
          if (i!=j) {
            const Particle& p_i = part[i];
            const Particle& p_j = part[j];
            double dphi = deltaPhi(p_i, p_j);
            double dX = Xj[i] - Xj[j];
            c_E   += cos(omega*dX - dphi);
          }
        }
      }
      return c_E/part.size();
    }


    void fillS(AIDA::IHistogram1D* h, const ParticleVector& part, double weight, std::vector<double> Xj, bool SE=true) {
      for (int i=0; i< h->axis().bins(); i++) {
        double x = h->axis().binLowerEdge(i);
        x += 0.5 * (h->axis().binUpperEdge(i) - h->axis().binLowerEdge(i));
        double y;
        if (SE) y = getSE(part, Xj, x);
        else    y = getSeta(part, x);
        h->fill(x, y*weight);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
         
      double weight = event.weight(); 
      // Charged fs
      const ChargedFinalState& cfs100 = applyProjection<ChargedFinalState>(event, "CFS100");
      const ParticleVector    part100 = cfs100.particlesByEta();
      
      const ChargedFinalState& cfs500 = applyProjection<ChargedFinalState>(event, "CFS500");
      const ParticleVector&   part500 = cfs500.particlesByEta();
      
      // The most first the pTmax < 10 and pT > 100 MeV part
    


      if (part100.size() > 10) {
        double ptmax100                 = cfs100.particlesByPt()[0].momentum().pT()/GeV;
        if (ptmax100 < 10) {
          std::vector<double> Xj100 = getXj(part100);
          fillS(_sE_10_100, part100, weight, Xj100, true);
          fillS(_sEta_10_100, part100, weight, Xj100, false);
        }
      }

      // Simple counter for the pTmax < 1.0 GeV case
      int nsmallpT = 0;
      ParticleVector smallpT;
      foreach (const Particle& p, cfs100.particlesByPt()) {
        if (p.momentum().pT()/GeV < 1.0) nsmallpT++;
        smallpT.push_back(p);
      }
      
      if (nsmallpT > 10) {
        std::vector<double> XjsmallpT = getXj(smallpT);
        fillS(_sE_1_100, smallpT, weight, XjsmallpT, true);
        fillS(_sEta_1_100, smallpT, weight, XjsmallpT, false);
      }
      
      if (part500.size() > 10) {
        double ptmax500                 = cfs500.particlesByPt()[0].momentum().pT()/GeV;
        if (ptmax500 < 10) {
          std::vector<double> Xj500 = getXj(part500);
          fillS(_sE_10_500, part500, weight, Xj500, true);
          fillS(_sEta_10_500, part500, weight, Xj500, false);
        }
      }
    }
      
    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sE_10_100, 1.0/(sumOfWeights()*_sE_10_100->axis().bins()));
      scale(_sE_1_100 , 1.0/(sumOfWeights()*_sE_1_100 ->axis().bins()));
      scale(_sE_10_500, 1.0/(sumOfWeights()*_sE_10_500->axis().bins()));  
                    
      scale(_sEta_10_100, 1.0/(sumOfWeights()*_sEta_10_100->axis().bins()));
      scale(_sEta_1_100 , 1.0/(sumOfWeights()*_sEta_1_100 ->axis().bins()));
      scale(_sEta_10_500, 1.0/(sumOfWeights()*_sEta_10_500->axis().bins()));
    } 

    //@}


  private:

    AIDA::IHistogram1D* _sE_10_100;
    AIDA::IHistogram1D* _sE_1_100;
    AIDA::IHistogram1D* _sE_10_500;  
                
    AIDA::IHistogram1D* _sEta_10_100;
    AIDA::IHistogram1D* _sEta_1_100;
    AIDA::IHistogram1D* _sEta_10_500;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1091481);

}
