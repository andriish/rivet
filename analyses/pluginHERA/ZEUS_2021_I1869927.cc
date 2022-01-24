// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#include "Rivet/AnalysisHandler.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

namespace Rivet {


  /// @brief dN/dNch, dN/dpT, dN/deta, and two-particle azimuthal correlations vs DeltaEta and MeanPt for charged hadrons in photoproduction
  class ZEUS_2021_I1869927 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ZEUS_2021_I1869927);


    /// @name Analysis methods
    ///@{
    const int MultLowLimit = 20;
   
    const vector<int> ParentDecayIDs = {2212, 22, 311, 11, 2112, 13, 130, 211, 321, 
	3322, 3122, 3312, 3112, 310, 3334, 3222};

    // Remove particles if they were decay products of the above particles
    bool UnwantedDecay( Particle p ) {
	for( auto el : ParentDecayIDs ) {
	    if( p.hasParent( el ) ) { return true; }
	}
	return false;
    }



    void init() {

      // Particles to be used for analysis
      // 0.1 < pT < 5.0 GeV, -1.5 < eta < 2.0
      // charged pions, kaons, protons, Xi, Sigma, Omega
	const ChargedFinalState cfs( 
		Cuts::pT > 0.1*GeV && Cuts::pT < 5.0*GeV && 
		Cuts::eta > -1.5 && Cuts::eta < 2.0 && 
		( Cuts::abspid == 211 || Cuts::abspid == 321 || Cuts::abspid == 2212 || 
		  Cuts::abspid == 3312 || Cuts::abspid == 3222 || Cuts::abspid == 3112 || Cuts::abspid == 3334 )
		);

      declare(cfs, "SelectedParticles");
      
      // Automatic binning from HEP data yoda file: ZEUS_2021_I1869927.yoda
      // HEP data file only has the bin centers defined so this method doesn't work
      // format: tableNumber, xAxisNumber, yAxisNumber
      // Only 1 type of X and Y axis per table so these are always 1, 1
      //book( _h_Nch, 9, 1, 1 );
      //book( _h_pT, 10, 1, 1 );
      //book( _h_eta, 11, 1, 1 );
      //book( _h_c12_dEta, 12, 1, 1 ); 
      //book( _h_c22_dEta, 13, 1, 1 ); 
      //book( _h_c12_mPt, 14, 1, 1 ); 
      //book( _h_c22_mPt, 15, 1, 1 ); 


      // Manually specified bins
      book( _h_Nch, "Nch", {19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5} );
      book( _h_pT, "pT", {0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5.0} );
      book( _h_eta, "eta", {-1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0} );
      book( _h_c12_dEta, "c12_dEta", {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.5} ); 
      book( _h_c22_dEta, "c22_dEta", {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.5} ); 
      book( _h_c12_mPt, "c12_mPt", {0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.5} ); 
      book( _h_c22_mPt, "c22_mPt", {0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.5} ); 

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
   
	// calculate total Nch
	int Nch = 0;
        const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event,"SelectedParticles");
        for(const Particle& p1 : cfs.particles()) {
	   if( UnwantedDecay( p1 ) ) { continue; }
           Nch++;
        }
	
	// Cut on Nch
	if( Nch < MultLowLimit ) { return; }

	_h_Nch->fill( Nch );

	//
	// Start main analysis
	//

	// 1st particle loop
	for(const Particle& p1 : cfs.particles()) {
	    if( UnwantedDecay( p1 ) ) { continue; }

	   _h_pT->fill(p1.pT()/GeV);
           _h_eta->fill(p1.eta());

	   // 2nd particle loop
	   for(const Particle& p2 : cfs.particles()) {
	       if( UnwantedDecay( p2 ) ) { continue; }

	       double dEta = fabs( p1.eta() - p2.eta() );
	       double mPt = ( p1.pT() - p2.pT() )/2.;

	       double c1 = cos( p1.phi() - p2.phi() );
	       double c2 = cos( 2*(p1.phi() - p2.phi()) );
	       
	       _h_c12_dEta->fill( dEta, c1 );
	       _h_c22_dEta->fill( dEta, c2 );
	       
	       _h_c12_mPt->fill( mPt, c1 );
	       _h_c22_mPt->fill( mPt, c2 );
	   }
 
	}

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize( _h_Nch ); // normalize to unity
      normalize( _h_eta ); // normalize to unity
      normalize( _h_pT ); // normalize to unity

      cout<<"Done"<<endl;
    }

    ///@}

      private:

    Histo1DPtr _h_Nch;
    Histo1DPtr _h_pT;
    Histo1DPtr _h_eta;
    
    Profile1DPtr _h_c12_dEta;
    Profile1DPtr _h_c22_dEta;
    Profile1DPtr _h_c12_mPt;
    Profile1DPtr _h_c22_mPt;
    

  };


  DECLARE_RIVET_PLUGIN(ZEUS_2021_I1869927);

}
