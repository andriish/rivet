// -*- C++ -*-
// -*- C++ -*-
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {


  /// @brief Identified particles in p--Pb @ 5 TeV
  class ALICE_2014_I1244523 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2014_I1244523);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      
      // The centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(),
           "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Define the cuts for the analysis:
      // pPb Collision has a centre of mass system shift of +0.465
      // They study -0.5 < yCoM < 0.0 -> -0.035 < y < 0.465
      const Cut& cut = Cuts::rap > -0.035 && Cuts::rap < 0.465;
      const ALICE::PrimaryParticles fs(cut);
      addProjection(fs,"FS");

      
      // The event trigger.
      declare(ALICE::V0AndTrigger(), "V0-AND");


      // Cut to figure out which multiplicity class the event lies in
      // Calculated via -0.5 < etaLab < 0.5
      const Cut& cutMultClass = Cuts::abseta<0.5;
      const ChargedFinalState cfsMultClass(cutMultClass);
      addProjection(cfsMultClass,"CFSMC");

      // Book histograms
      // Multiplicity classes for the analysis
      // Store the histograms in an array for easy event loop analysis
      // 2 sets of d0x for each of these histograms
      for (int ihist = 0; ihist < 4; ++ihist) {
        _histPipT[ihist]           = bookHisto1D(1,1,1+ihist);
        _histKpT[ihist]            = bookHisto1D(3,1,1+ihist);
        _histKS0pT[ihist]          = bookHisto1D(5,1,1+ihist);
        _histProtonpT[ihist]       = bookHisto1D(7,1,1+ihist);
        _histLambdapT[ihist]       = bookHisto1D(9,1,1+ihist);

        _histNPi4K[ihist]          = bookHisto1D("TMP/NPi4K",refData(11,1,1+ihist));
        _histNPi4Proton[ihist]     = bookHisto1D("TMP/NPi4Pro",refData(13,1,1+ihist));
        _histNKS04Lambda[ihist]    = bookHisto1D("TMP/NK204Lam",refData(15,1,1+ihist));

        _histNK[ihist]             = bookHisto1D("TMP/NK",refData(11,1,1+ihist));
        _histNProton[ihist]        = bookHisto1D("TMP/NPro",refData(13,1,1+ihist));
        _histNLambda[ihist]        = bookHisto1D("TMP/NLam",refData(15,1,1+ihist));

        _histNKtoPi[ihist]         = bookScatter2D(11,1,1+ihist);
        _histNProtontoPi[ihist]    = bookScatter2D(13,1,1+ihist);
        _histNLambdatoKS0[ihist]   = bookScatter2D(15,1,1+ihist);

      }
      for (int i = 0; i < 3; ++i) {
        int ihist = i + 4;
        _histPipT[ihist]          = bookHisto1D(2,1,1+i);
        _histKpT[ihist]           = bookHisto1D(4,1,1+i);
        _histKS0pT[ihist]         = bookHisto1D(6,1,1+i);
        _histProtonpT[ihist]      = bookHisto1D(8,1,1+i);
        _histLambdapT[ihist]      = bookHisto1D(10,1,1+i);

        _histNPi4K[ihist]         = bookHisto1D("TMP/NPi4K",refData(12,1,1+i));
        _histNPi4Proton[ihist]    = bookHisto1D("TMP/NPi4Pro",refData(14,1,1+i));
        _histNKS04Lambda[ihist]   = bookHisto1D("TMP/NK204Lam",refData(16,1,1+i));

        _histNK[ihist]            = bookHisto1D("TMP/NK",refData(12,1,1+i));
        _histNProton[ihist]       = bookHisto1D("TMP/NPro",refData(14,1,1+i));
        _histNLambda[ihist]       = bookHisto1D("TMP/NLam",refData(16,1,1+i));

        _histNKtoPi[ihist]        = bookScatter2D(12,1,1+i);
        _histNProtontoPi[ihist]   = bookScatter2D(14,1,1+i);
        _histNLambdatoKS0[ihist]  = bookScatter2D(16,1,1+i);
      }

      _histLambdaMeanpT           = bookProfile1D(17, 1, 1);
      _histProtonMeanpT           = bookProfile1D(18, 1, 1);
      _histKS0MeanpT              = bookProfile1D(19, 1, 1);
      _histKMeanpT                = bookProfile1D(20, 1, 1);
      _histPiMeanpT               = bookProfile1D(21, 1, 1);

      _histKtoPiYield             = bookScatter2D(22,1,1);
      _histProtontoPiYield        = bookScatter2D(22,1,2);
      _histLambdatoKS0Yield       = bookScatter2D(22,1,3);

      _histKYield                 = bookHisto1D("TMP/KY", refData(22,1,1));
      _histProtonYield            = bookHisto1D("TMP/PrY",refData(22,1,2));
      _histLambdaYield            = bookHisto1D("TMP/LY", refData(22,1,3));
      _histPiYield                = bookHisto1D("TMP/PiY",refData(22,1,1));
      _histKS0Yield               = bookHisto1D("TMP/KSY",refData(22,1,3));


      _multBinValues[0] = 45.0;
      _multBinValues[1] = 36.2;
      _multBinValues[2] = 30.5;
      _multBinValues[3] = 23.3;
      _multBinValues[4] = 16.1;
      _multBinValues[5] = 9.8;
      _multBinValues[6] = 4.3;

      _multBinCent[0]   = 2.5;
      _multBinCent[1]   = 7.5;
      _multBinCent[2]   = 15.0;
      _multBinCent[3]   = 30.0;
      _multBinCent[4]   = 50.0;
      _multBinCent[5]   = 70.0;
      _multBinCent[6]   = 90.0;


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Count the number of particles in an event
      // to choose the multiplicity class
      const double weight = event.weight();
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event,"CFSMC");
      
      // Event trigger.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")() ) vetoEvent;

      //Centrality
      const CentralityProjection& cent = apply<CentralityProjection>(event,"V0M");
      double c = cent();
      
      int centralityclass = -1;
      if(c > 0. && c <= 10) centralityclass = 0;
      if(c > 10. && c <= 20) centralityclass = 1;
      if(c > 20. && c <= 40) centralityclass = 2;
      if(c > 40. && c <= 60) centralityclass = 3;
      if(c > 60. && c <= 80) centralityclass = 4;
      if (centralityclass==-1) vetoEvent;
      
      const ALICE::PrimaryParticles& fs =
        apply<ALICE::PrimaryParticles>(event,"FS");
      foreach(const Particle& p, fs.particles()) {
        // Protections against MC generators decaying long lived particles
          switch (abs(p.pid())){
            case 211: // Pions for the ratios
              _histPipT[centralityclass]->fill(p.pT()/GeV, weight/(2.*M_PI*p.pT()/GeV));
              _histNPi4K[centralityclass]->fill(p.pT()/GeV, weight);
              _histNPi4Proton[centralityclass]->fill(p.pT()/GeV, weight);
              _histPiYield->fill(_multBinCent[centralityclass], weight);
              break;
            case 3312: // K +/-
              _histKpT[centralityclass]->fill(p.pT()/GeV, weight/(2.*M_PI*p.pT()/GeV));
              _histNK[centralityclass]->fill(p.pT()/GeV, weight);
              _histKYield->fill(_multBinCent[centralityclass], weight);
              break;
            case 3334: // Ks0
              _histKS0pT[centralityclass]->fill(p.pT()/GeV, 2.0*weight/(2.*M_PI*p.pT()/GeV));
              _histNKS04Lambda[centralityclass]->fill(p.pT()/GeV, weight);
              _histKS0Yield->fill(_multBinCent[centralityclass], weight);
              break;
            case 2212: // proton
              _histProtonpT[centralityclass]->fill(p.pT()/GeV, weight/(2.*M_PI*p.pT()/GeV));
              _histNProton[centralityclass]->fill(p.pT()/GeV, weight);
              _histProtonYield->fill(_multBinCent[centralityclass], weight);
              break;
            case 3122: // Lambda
              _histLambdapT[centralityclass]->fill(p.pT()/GeV, weight/(2.*M_PI*p.pT()/GeV));
              _histNLambda[centralityclass]->fill(p.pT()/GeV,weight);
              _histLambdaYield->fill(_multBinCent[centralityclass], weight);
              break;
          } // Particle Switch
      } // Particle loop
    }// End of init

    /// Normalise histograms etc., after the run
    void finalize() {

      // Do mean pT spectrum first, since the next loop modifies the histograms
      for (int i = 0; i < 7; i++){
        _histKMeanpT->fill( _multBinValues[i],     _histKpT[i]->xMean());
        _histKS0MeanpT->fill( _multBinValues[i],   _histKS0pT[i]->xMean());
        _histProtonMeanpT->fill(_multBinValues[i], _histProtonpT[i]->xMean());
        _histLambdaMeanpT->fill(_multBinValues[i], _histLambdapT[i]->xMean());
        _histPiMeanpT->fill(_multBinValues[i],     _histPipT[i]->xMean());
      }
      // Now we can scale the pT spectrum histograms
      for (int i = 0; i < 7; i++){

        divide(_histNK[i],      _histNPi4K[i],       _histNKtoPi[i]);
        divide(_histNProton[i], _histNPi4Proton[i],  _histNProtontoPi[i]);
        divide(_histNLambda[i], _histNKS04Lambda[i], _histNLambdatoKS0[i]);

        scale( _histNK[i],         1.0/sumOfWeights());
        scale( _histNProton[i],    1.0/sumOfWeights());
        scale( _histNLambda[i],    1.0/sumOfWeights());

        scale( _histPipT[i],       1.0/sumOfWeights());
        scale( _histKpT[i],        1.0/sumOfWeights());
        scale( _histKS0pT[i],      1.0/sumOfWeights());
        scale( _histProtonpT[i],   1.0/sumOfWeights());
        scale( _histLambdapT[i],   1.0/sumOfWeights());
      }

      divide(_histKYield,      _histPiYield,  _histKtoPiYield);
      divide(_histProtonYield, _histPiYield,  _histProtontoPiYield);
      divide(_histLambdaYield, _histKS0Yield, _histLambdatoKS0Yield);

    }

    /*
    int _calculateMultClass(const ChargedFinalState& cfs){

    }
    */

    //@}

private:
    // pT spectrum vs number produced (separated by multiplicity classes)
    Histo1DPtr         _histPipT[7];
    Histo1DPtr         _histKpT[7];
    Histo1DPtr         _histKS0pT[7];
    Histo1DPtr         _histProtonpT[7];
    Histo1DPtr         _histLambdapT[7];
    // pT vs ratio of yields (separated by multiplicity classes)
    Histo1DPtr         _histNPi4K[7];
    Histo1DPtr         _histNPi4Proton[7];
    Histo1DPtr         _histNKS04Lambda[7];
    Histo1DPtr         _histNK[7];
    Histo1DPtr         _histNProton[7];
    Histo1DPtr         _histNLambda[7];

    Scatter2DPtr       _histNKtoPi[7];
    Scatter2DPtr       _histNProtontoPi[7];
    Scatter2DPtr       _histNLambdatoKS0[7];
    // Multiplicity class vs mean pT
    // (1 graph for each species, K, KS0, p, lambda)
    Profile1DPtr       _histKMeanpT;
    Profile1DPtr       _histKS0MeanpT;
    Profile1DPtr       _histProtonMeanpT;
    Profile1DPtr       _histLambdaMeanpT;
    Profile1DPtr       _histPiMeanpT;

    Histo1DPtr         _histKYield;
    Histo1DPtr         _histProtonYield;
    Histo1DPtr         _histLambdaYield;
    Histo1DPtr         _histPiYield;
    Histo1DPtr         _histKS0Yield;

    Scatter2DPtr       _histKtoPiYield;
    Scatter2DPtr       _histProtontoPiYield;
    Scatter2DPtr       _histLambdatoKS0Yield;

    double             _multBinValues[7];
    double             _multBinCent[7];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2014_I1244523);


}
