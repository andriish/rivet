// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Correlators.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "YODA/Scatter2D.h"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_HM_FLOW : public CumulantAnalysis {
  public:

    /// Constructor
    ALICE_HM_FLOW() : CumulantAnalysis("ALICE_HM_FLOW"){};


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
 // Initialise and register projections
      // Declare the trigger projection.
      declare<ALICE::V0AndTrigger>(ALICE::V0AndTrigger(),"V0-AND");
      declare<ALICE::V0MMultiplicity>(ALICE::V0MMultiplicity(),"V0M");
      
      // The full central charged final state.
      const ChargedFinalState& cfs = ChargedFinalState(Cuts::abseta < 0.8 &&
        Cuts::pT > 0.2*GeV && Cuts::pT < 5.0*GeV);
      declare(cfs, "CFS");

      // The positive eta side used for rapidity gap.
      const ChargedFinalState& cfsp = ChargedFinalState(Cuts::eta > 0.5 && 
        Cuts::eta < 0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 3.0*GeV);
      declare(cfsp, "CFSP");
      // ..negative ditto.
      const ChargedFinalState& cfsn = ChargedFinalState(Cuts::eta < -0.5 && 
        Cuts::eta > -0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 3.0*GeV);
      declare(cfsn, "CFSN");

      // We need a larger gap for v2.
      const ChargedFinalState& cfspl = ChargedFinalState(Cuts::eta > 0.7 && 
        Cuts::eta < 0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 3.0*GeV);
      declare(cfsp, "CFSPL");
      // ..negative ditto.
      const ChargedFinalState& cfsnl = ChargedFinalState(Cuts::eta < -0.7 && 
        Cuts::eta > -0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 3.0*GeV);
      declare(cfsn, "CFSNL");

      h_v22 = bookScatter2D("v22",120,0,120);
      h_v32 = bookScatter2D("v32",120,0,120);
      h_v42 = bookScatter2D("v42",120,0,120);

      h_v22gap = bookScatter2D("v22gap",120,0,120);
      h_v32gap = bookScatter2D("v32gap",120,0,120);
      h_v42gap = bookScatter2D("v42gap",120,0,120);

      vector<double> pTbins = {0.2,0.5,1.0,2.0,3.0,4.0,5.0,7.0,9.0,11.0,15.0};

      h_v22gappT = bookScatter2D("v22gappT",pTbins,"v22gappT","pT","v2");
      h_v32gappT = bookScatter2D("v32gappT",pTbins,"v32gappT","pT","v3");
      h_v42gappT = bookScatter2D("v42gappT",pTbins,"v42gappT","pT","v4");

      h_v24 = bookScatter2D("v24",120,0,120);
      
      h_c22 = bookScatter2D("c22",120,0,120);
      h_c32 = bookScatter2D("c32",120,0,120);
      h_c42 = bookScatter2D("c42",120,0,120);

      h_c22gap = bookScatter2D("c22gap",120,0,120);
      h_c32gap = bookScatter2D("c32gap",120,0,120);
      h_c42gap = bookScatter2D("c42gap",120,0,120);

      h_c24 = bookScatter2D("c24",120,0,120);

      h_ec22 = bookScatter2D("ec22",120,0,120);      
      h_ec24 = bookScatter2D("ec24",120,0,120);      

      h_ec22fm = bookScatter2D("ec22fm",120,0,120);      
      h_ec24fm = bookScatter2D("ec24fm",120,0,120);      
      
      h_c22fm = bookScatter2D("c22fm",120,0,120);      
      h_c24fm = bookScatter2D("c24fm",120,0,120);      
      
      // Corresponding event averaged correlators.
      // Integrated, no gap.
      ec22 = bookECorrelator<2,2>("ec22",h_v22);
      ec32 = bookECorrelator<3,2>("ec32",h_v32);
      ec42 = bookECorrelator<4,2>("ec42",h_v42);
      ec24 = bookECorrelator<2,4>("ec24",h_v24);

      ec22fm = bookECorrelator<2,2>("ec22fm",h_c22fm);
      ec24fm = bookECorrelator<2,2>("ec24fm",h_c24fm);

      // ... with gap.
      ec22gap = bookECorrelatorGap<2,2>("ec22gap",h_v22gap);
      ec32gap = bookECorrelatorGap<3,2>("ec32gap",h_v32gap);
      ec42gap = bookECorrelatorGap<4,2>("ec42gap",h_v42gap);
      // ... pT binned
      ec22gappT = bookECorrelatorGap<2,2>("ec22gappT",h_v22gappT);
      ec32gappT = bookECorrelatorGap<3,2>("ec32gappT",h_v32gappT);
      ec42gappT = bookECorrelatorGap<4,2>("ec42gappT",h_v42gappT); 
   

      // Maximal N and P.
      pair<int, int> max = getMaxValues(); 
      // Declare correlator projections.
      declare(Correlators(cfs, max.first, max.second),"Correlators");
      declare(Correlators(cfsp, max.first, max.second, pTbins),"CorrelatorsPos");
      declare(Correlators(cfsn, max.first, max.second, pTbins),"CorrelatorsNeg");
      declare(Correlators(cfspl, max.first, max.second),"CorrelatorsPosL");
      declare(Correlators(cfsnl, max.first, max.second),"CorrelatorsNegL");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Event trigger.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")() ) vetoEvent;
      
      // High multiplicity trigger.
      if (apply<ALICE::V0MMultiplicity>(event, "V0M")() < 160.) vetoEvent;
      

      const double w = event.weight();
      
      // Projections.
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event,"CFS");

      // The cumulants projections.
      const Correlators& c = applyProjection<Correlators>(event,"Correlators");
      const Correlators& cn = applyProjection<Correlators>(event,"CorrelatorsPos");
      const Correlators& cp = applyProjection<Correlators>(event,"CorrelatorsNeg");
      const Correlators& cnl = applyProjection<Correlators>(event,"CorrelatorsPosL");
      const Correlators& cpl = applyProjection<Correlators>(event,"CorrelatorsNegL");

      // Number of charged particles.
      const double nch = cfs.particles().size();

      const double nchFwd = apply<ALICE::V0MMultiplicity>(event,"V0M")();

      ec22->fill(nch, c, w);
      ec32->fill(nch, c, w);
      ec42->fill(nch, c, w);
      ec24->fill(nch, c, w);

      ec22fm->fill(nchFwd, c, w);
      ec24fm->fill(nchFwd, c, w);

      ec22gap->fill(nch, cpl, cnl, w);
      ec32gap->fill(nch, cp, cn, w);
      ec42gap->fill(nch, cp, cn, w);
      
      ec22gappT->fill(cp, cn, w);
      ec32gappT->fill(cp, cn, w);
      ec42gappT->fill(cp, cn, w);
    }



    /// Normalise histograms etc., after the run
    void finalize() {

      cnTwoInt(h_c22,ec22);
      cnTwoInt(h_c32,ec32);
      cnTwoInt(h_c42,ec42);

      cnTwoInt(h_c22gap,ec22gap);
      cnTwoInt(h_c32gap,ec32gap);
      cnTwoInt(h_c42gap,ec42gap);
      
      cnFourInt(h_c24, ec22, ec24);
      
      vnTwoInt(h_v22,ec22);
      vnTwoInt(h_v32,ec32);
      vnTwoInt(h_v42,ec42);

      vnTwoInt(h_v22gap,ec22gap);
      vnTwoInt(h_v32gap,ec32gap);
      vnTwoInt(h_v42gap,ec42gap);
      
      vnFourInt(h_v24, ec22, ec24);

      vnTwoDiff(h_v22gappT, ec22gappT);
      vnTwoDiff(h_v32gappT, ec32gappT);
      vnTwoDiff(h_v42gappT, ec42gappT);

      corrPlot(h_ec22, ec22);
      corrPlot(h_ec24, ec24);

      cnTwoInt(h_c22fm, ec22fm);
      cnFourInt(h_c24fm, ec22fm, ec24fm);
      corrPlot(h_ec22fm, ec22fm);
      corrPlot(h_ec24fm, ec24fm);



    }
    //@}


    /// @name Histograms
    //@{
      Scatter2DPtr h_v22;
      Scatter2DPtr h_v32;
      Scatter2DPtr h_v42;

      Scatter2DPtr h_v22gap;
      Scatter2DPtr h_v32gap;
      Scatter2DPtr h_v42gap;
      
      Scatter2DPtr h_v22gappT;
      Scatter2DPtr h_v32gappT;
      Scatter2DPtr h_v42gappT;
      
      Scatter2DPtr h_v24;
      
      Scatter2DPtr h_c22;
      Scatter2DPtr h_c32;
      Scatter2DPtr h_c42;
      
      Scatter2DPtr h_c22fm;
      Scatter2DPtr h_c24fm;

      Scatter2DPtr h_c22gap;
      Scatter2DPtr h_c32gap;
      Scatter2DPtr h_c42gap;

      Scatter2DPtr h_c24;
      
      // Test histos
      Scatter2DPtr h_ec22;
      Scatter2DPtr h_ec24;    
      Scatter2DPtr h_ec22fm;
      Scatter2DPtr h_ec24fm;

      // Event correlators.
      ECorrPtr ec22;
      ECorrPtr ec32;
      ECorrPtr ec42;
      // ... with gap.
      ECorrPtr ec22gap;
      ECorrPtr ec32gap;
      ECorrPtr ec42gap;

      ECorrPtr ec22gappT;
      ECorrPtr ec32gappT;
      ECorrPtr ec42gappT;
      // Four particle.
      ECorrPtr ec24;

      ECorrPtr ec22fm;
      ECorrPtr ec24fm;




    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_HM_FLOW);


}
