// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet{

  class ATLAS_2010_S8817804;
  
  AnalysisBuilder<ATLAS_2010_S8817804> plugin_ATLAS_2010_S8817804;
  
  /// @brief ATLAS inclusive jet pT spectrum, di-jet Mass and di-jet chi
  
  class ATLAS_2010_S8817804: public Analysis{
    
  public:
    
    ATLAS_2010_S8817804(): Analysis("ATLAS_2010_S8817804"), _rapMax(2.8), _leadPTCut(60. * GeV), _asymmetryCut(log(30.)), _ySumCut(2.2){
      setBeams(PROTON, PROTON);
      setNeedsCrossSection(true);
    }
    
    void init(){
      FinalState fs;
      addProjection(fs, "FinalState");
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.6), "AntiKT06");
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.4), "AntiKT04");      

      double ybins[] = {0.0, 0.3, 0.8, 1.2, 2.1, 2.8};

      double massBinsForChi[] = { 340. * GeV, 520. * GeV, 800. * GeV, 1200. * GeV};//, 1700. * GeV, 2500. * GeV};
      
      int ptDsOffset = 1;
      int massDsOffset = 11;
      int chiDsOffset = 21;
      
      for(Jet_Alg alg = AKT4; alg != ALG_END; alg = Jet_Alg(alg + 1)){
        for(int highbin = 0, lowbin = highbin++; highbin != sizeof(ybins) / sizeof(double); ++lowbin, ++highbin){

          _pTHistos[alg].addHistogram(ybins[lowbin], ybins[highbin],
                                      bookHistogram1D(lowbin + ptDsOffset, 1, 1));
        }
        
        ptDsOffset = 6;
        
        for(int highbin = 0, lowbin = highbin++; highbin != sizeof(ybins) / sizeof(double); ++lowbin, ++highbin){

          _massVsY[alg].addHistogram(ybins[lowbin], ybins[highbin],
                                     bookHistogram1D(massDsOffset + lowbin, 1, 1));
        }
       
        massDsOffset = 16;
        
        for(int highbin=0, lowbin = highbin++; highbin != sizeof(massBinsForChi) / sizeof(double); ++lowbin, ++highbin){

          _chiVsMass[alg].addHistogram(massBinsForChi[lowbin], massBinsForChi[highbin],
                                       bookHistogram1D(chiDsOffset + lowbin, 1, 1));
          
        }
        
        chiDsOffset = 24;
        
      }
            
      return;
    }
    
    void analyze(const Event &evt){
        
      const Jets &kt6Jets = applyProjection<FastJets>(evt, "AntiKT06").jets(30.*GeV);
      const Jets &kt4Jets = applyProjection<FastJets>(evt, "AntiKT04").jets(30.*GeV);
      
      const Jets* jetAr[2];
      jetAr[AKT4] = &kt4Jets;
      jetAr[AKT6] = &kt6Jets;
                  
      for(Jet_Alg alg = AKT4; alg != ALG_END; alg = Jet_Alg(alg + 1)){
        
        const Jet *jet1 = 0, *jet2 = 0;
        double pt1 = -1.;
        double pt2 = -1.;
        
        for(Jets::const_iterator jet = jetAr[alg]->begin();
            jet != jetAr[alg]->end(); ++jet){
          double pT = jet->momentum().pT();
          double fabsy = fabs(jet->momentum().rapidity());
          _pTHistos[alg].fill(fabsy, pT, evt.weight());
          
          if(fabsy < _rapMax){
            if(pT > pt1){
              if(jet1 != 0){
                jet2 = jet1;
                pt2 = jet2->momentum().pT();
              }
              jet1 = &(*jet);
              pt1 = pT;
            }else if(pT > pt2){
              jet2 = &(*jet);
              pt2 = pT;
            }
          }
        }
        
        if( jet1 == 0 || jet2 == 0 || pt1 < _leadPTCut){
          continue;
        }
        
        FourMomentum mom1 = jet1->momentum();
        FourMomentum mom2 = jet2->momentum();
        FourMomentum total = mom1 + mom2;
        double rap1 = mom1.rapidity();
        
        double absRap1 = fabs(rap1);
        
        double rap2 = mom2.rapidity();
        double absRap2 = fabs(rap2);
        double ymax = (absRap1 > absRap2)? absRap1: absRap2;
        double chi = exp(fabs(rap1 - rap2));
        double mass = total.mass();
        
        if(fabs(rap1 + rap2) < _ySumCut && fabs(rap1 - rap2) < _asymmetryCut){
          _chiVsMass[alg].fill(mass, chi, evt.weight());
        }
        
        _massVsY[alg].fill(ymax, mass, evt.weight());
        
      }
      
      return;
    }
    
    void finalize(){
            
      const double xs = crossSectionPerEvent() / picobarn;
            
      for(Jet_Alg alg = AKT4; alg != ALG_END; alg = Jet_Alg(alg + 1)){
        _pTHistos[alg].scale(xs, this);
        _massVsY[alg].scale(xs, this);
        _chiVsMass[alg].scale(xs, this);
      }
  
      return;
    }
    
  private:
    
    // The max rapidity of a jet used in this analysis
    double _rapMax;
    // The pT cut on the leading jet 
    double _leadPTCut;
    // maximum rapidity asymmetry between two jets in the chi plots
    double _asymmetryCut;
    // maximum sum of rapidities of two jets in the chi plot
    double _ySumCut;
    
    enum Jet_Alg{AKT4=0, AKT6=1, ALG_END};
        
    // The inclusive pT spectrum for akt6 and akt4 jets (array index is jet type from enum above)
    BinnedHistogram<double> _pTHistos[2];
    
    // The di-jet mass spectrum binned in rapidity for akt6 and akt4 jets (array index is jet type from enum above)
    BinnedHistogram<double> _massVsY[2];
    
    // The di-jet chi distribution binned in mass for akt6 and akt4 jets (array index is jet type from enum above)
    BinnedHistogram<double> _chiVsMass[2];
    
    
  };
}
