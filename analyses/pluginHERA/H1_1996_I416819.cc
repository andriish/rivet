// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class H1_1996_I416819 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(H1_1996_I416819);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
		
	
       Histo1DPtr dummy;
	 _h_f2.add( 1.0,     2.3,book(dummy,1,1,1));
	 //_h_f2.add( 2.2,     2.85,book(dummy,"TMP/dummy1", refData(1,1,1))); // remove this region for better coverage   
	 _h_f2.add( 2.85,    4.3,book(dummy,4,1,1));

	 _h_f2.add( 4.3,     5.9,book(dummy,5,1,1));
	 _h_f2.add( 5.9,     7.2,book(dummy,6,1,1));

	 _h_f2.add( 7.2,    10.1,book(dummy,7,1,1));


	 // _h_f2.add( 10.1,   11.5,book(dummy,"TMP/dummy2", refData(7,1,1)));   // remove this region for better coverage   
	 _h_f2.add( 11.5,   12.5,book(dummy,8,1,1));
	 _h_f2.add( 12.5,   18.3,book(dummy,9,1,1));
	 _h_f2.add( 18.3,   22., book(dummy,10,1,1));
	 _h_f2.add( 22.,    29,  book(dummy,11,1,1));
	 _h_f2.add( 29.,    42., book(dummy,12,1,1));
	 _h_f2.add( 42.,    50., book(dummy,13,1,1));
	 _h_f2.add( 50.,    73., book(dummy,14,1,1));
	 _h_f2.add( 73.,    115.,book(dummy,15,1,1));
	 _h_f2.add( 115.,   130.,book(dummy,16,1,1));
	 _h_f2.add( 130.,   180.,book(dummy,17,1,1));
	 _h_f2.add( 180.,   233.,book(dummy,18,1,1));
	 _h_f2.add( 233.,   280.,book(dummy,19,1,1));
	 _h_f2.add( 280.,   455.,book(dummy,20,1,1));
	 _h_f2.add( 455.,   558.,book(dummy,21,1,1));
	 _h_f2.add( 558.,   780.,book(dummy,22,1,1));
	 _h_f2.add( 780.,   830.,book(dummy,23,1,1));
	 _h_f2.add( 830.,  1290.,book(dummy,24,1,1));
	 _h_f2.add( 1290., 4000.,book(dummy,25,1,1));
	 _h_f2.add( 4000., 6700.,book(dummy,26,1,1));


/*
        book(_hist_Q2_10, "Q2_10",100,1., 11.0);
        book(_hist_Q2_100, "Q2_100",100,10., 100.0);
        book(_hist_Q2_1000, "Q2_1000",100,100., 1000.0);
        book(_hist_Q2_2000, "Q2_2000",100,800., 5000.0);
        book(_hist_Q2_3000, "Q2_3000",100,3000., 10000.0);
*/

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
      //const DISLepton& dl = applyProjection<DISLepton>(event,"Lepton");

      // Get the DIS kinematics
      double x  = dk.x();
      double y = dk.y();
      double Q2 = dk.Q2()/GeV;
	
	// Flux factor
	const double alpha = 7.29927e-3;
	double F = x*pow(Q2,2.)/(2.0*M_PI*pow(alpha,2.)*(1.0+pow((1.-y),2.)));
/*
	_hist_Q2_10-> fill(Q2) ;
	_hist_Q2_100-> fill(Q2) ;
	_hist_Q2_1000-> fill(Q2) ;
	_hist_Q2_2000-> fill(Q2) ;
	_hist_Q2_3000-> fill(Q2) ;
*/
	_h_f2.fill(Q2,x,F); // wypelniamy histogram x w skali Q2

    }


    /// Normalise histograms etc., after the run
    void finalize() {
	double gev2nb =0.389e6;
      double scalefactor=crossSection()/nanobarn/sumOfWeights()/gev2nb ;
      // with _h_f2.scale also q2 bin width is scaled
      _h_f2.scale(scalefactor, this);

    }

    ///@}


    /// @name Histograms
    ///@{
	BinnedHistogram _h_f2;
      Histo1DPtr _hist_Q2_10,_hist_Q2_100,_hist_Q2_1000,_hist_Q2_2000,_hist_Q2_3000;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(H1_1996_I416819);

}
