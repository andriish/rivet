// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "LWH/Histogram1D.h"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Particle.hh"

//#include <iostream>

namespace Rivet {


  class ATLAS_2012_I1188891 : public Analysis {

  public:


    ATLAS_2012_I1188891()
      : Analysis("ATLAS_2012_I1188891")
    {
    }



  public:

    void init() {

      const FinalState fs;
      FastJets fj04(fs,  FastJets::ANTIKT, 0.4);
      addProjection(fj04, "AntiKT04");

      //double xbins[7] = { 40.,60.,80.,120.,160.,250.,500. };
      vector<double> xbins;
      xbins.push_back(40.);
      xbins.push_back(60.);
      xbins.push_back(80.);
      xbins.push_back(120.);
      xbins.push_back(160.);
      xbins.push_back(250.);
      xbins.push_back(500.);

      string histotitle[7]={"BBfraction","BCfraction","CCfraction","BUfraction","CUfraction","UUfraction","Total"};
      for (int i = 0 ; i < 7 ; i++){
        _fracHistos[i] = bookHistogram1D(histotitle[i],xbins);
      }
    }


    void analyze(const Event& event) {

      double weight    = event.weight();
      double weight100 = event.weight() * 100.;  //to get results in %

      //keeps jets with pt>20 geV and ordered in decreasing pt
      Jets jetAr = applyProjection<FastJets>(event, "AntiKT04").jetsByPt(20*GeV);

      int flav[2]={1,1};
      vector<FourMomentum> leadjets;

      //get b/c-hadrons
      std::vector<HepMC::GenParticle*> B_hadrons;
      std::vector<HepMC::GenParticle*> C_hadrons;
      std::vector<HepMC::GenParticle*> allParticles = particles(event.genEvent());
      for(unsigned int i = 0; i < allParticles.size(); i++) {
        GenParticle* p = allParticles.at(i);
        if(p->momentum().perp()*GeV < 5) continue;
        if ( (Rivet::PID::isHadron ( p->pdg_id() ) &&
              Rivet::PID::hasBottom( p->pdg_id() )    ) ) {
          B_hadrons.push_back(p);
        }
        if ( (Rivet::PID::isHadron( p->pdg_id() ) &&
              Rivet::PID::hasCharm( p->pdg_id() )    ) ) {
          C_hadrons.push_back(p);
        }
      }

      //select dijet
      foreach (const Jet& jet, jetAr) {

        const double pT   = jet.momentum().pT();
        const double absy = fabs(jet.momentum().rapidity());

        bool isBjet = false;
        //not using this
        //isBjet = jet.containsBottom();
        foreach(HepMC::GenParticle* b, B_hadrons) {
          FourMomentum hadron = b->momentum();
          double hadron_jet_dR = deltaR(jet.momentum(), hadron);
          if(hadron_jet_dR < 0.3) isBjet = true;
        }

        bool isCjet = false;
        //bool isCjet = jet.containsCharm();
        foreach(HepMC::GenParticle* c, C_hadrons) {
          FourMomentum hadron = c->momentum();
          double hadron_jet_dR = deltaR(jet.momentum(), hadron);
          if(hadron_jet_dR < 0.3) isCjet = true;
        }

        int jetflav=1;
        if      (isBjet)jetflav=5;
        else if (isCjet)jetflav=4;

        if (absy <= 2.1 && leadjets.size() < 2) {

          if (pT > 500*GeV) continue;
          if ((leadjets.empty() && pT < 40*GeV) || pT < 20*GeV) continue;

          leadjets.push_back(jet.momentum());

          if (leadjets.size()==1) flav[0] = jetflav;
          if (leadjets.size()==2) flav[1] = jetflav;
        }
      }

      if (leadjets.size() < 2) vetoEvent;


      double pBinsLJ[7] = {40.,60.,80.,120.,160.,250.,500.};
      int    iPBinLJ = -1;

      for (int k = 0 ; k < 7 ; k++) {
        if (leadjets[0].pT() > pBinsLJ[k]*GeV) iPBinLJ=k;
        else break;
      }

      bool c_ljpt  = (iPBinLJ != -1);
      bool c_nljpt = leadjets[1].pT() > 20*GeV;
      bool c_dphi  = fabs( deltaPhi(leadjets[0],leadjets[1]) ) > 2.1;
      bool isDijet = c_ljpt & c_nljpt & c_dphi;
      if (!isDijet) vetoEvent;

      _fracHistos[6]->fill(leadjets[0].pT(), weight);

      if (flav[0]==5 && flav[1]==5)                                  // BB dijet
        _fracHistos[0]->fill(leadjets[0].pT(), weight100);

      if ((flav[0]==5 && flav[1]==4) || (flav[0]==4 && flav[1]==5))  // BC dijet
        _fracHistos[1]->fill(leadjets[0].pT(), weight100);

      if (flav[0]==4 && flav[1]==4)                                  // CC dijet
        _fracHistos[2]->fill(leadjets[0].pT(), weight100);

      if ((flav[0]==5 && flav[1]==1) || (flav[0]==1 && flav[1]==5))  // B-light dijet
        _fracHistos[3]->fill(leadjets[0].pT(), weight100);

      if ((flav[0]==4 && flav[1]==1) || (flav[0]==1 && flav[1]==4))  // C-light dijet
        _fracHistos[4]->fill(leadjets[0].pT(), weight100);

      if (flav[0]==1 && flav[1]==1)                                  // light-light dijet
        _fracHistos[5]->fill(leadjets[0].pT(), weight100);
    }


    void finalize() {
      AIDA::IHistogramFactory& hf = histogramFactory();

      hf.divide( histoPath(1,1,1), *_fracHistos[0], *_fracHistos[6]);
      hf.divide( histoPath(2,1,1), *_fracHistos[1], *_fracHistos[6]);
      hf.divide( histoPath(3,1,1), *_fracHistos[2], *_fracHistos[6]);
      hf.divide( histoPath(4,1,1), *_fracHistos[3], *_fracHistos[6]);
      hf.divide( histoPath(5,1,1), *_fracHistos[4], *_fracHistos[6]);
      hf.divide( histoPath(6,1,1), *_fracHistos[5], *_fracHistos[6]);

      // remove temporary histograms
      for (size_t i=0; i<7; i++) {
        hf.destroy(_fracHistos[i]);
      }
    }

  private:

    AIDA::IHistogram1D * _fracHistos[7];
  };

  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1188891);

}
