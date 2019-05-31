// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_D_Dalitz : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_D_Dalitz);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _h_plus1    = bookHisto1D("h_plus1"   ,200,0.,3. );
      _h_minus1   = bookHisto1D("h_minus1"  ,200,0.,3.2 );
      _h_pipi1    = bookHisto1D("h_pipi1"   ,200,0.,2. );
      _h_minus2   = bookHisto1D("h_minus2"  ,200,0.,3.2 );
      _h_neutral2 = bookHisto1D("h_neutral2",200,0.,3.2 );
      _h_pipi2    = bookHisto1D("h_pipi2"   ,200,0.,2. );
      _h_Kpilow3   = bookHisto1D("h_Kpilow3"  ,200,0.,2. );
      _h_Kpihigh3  = bookHisto1D("h_Kpihigh3" ,200,0.,3.2 );
      _h_Kpiall3   = bookHisto1D("h_Kpiall3"  ,200,0.,3. );
      _h_pipi3     = bookHisto1D("h_pipi3"    ,200,0.,2. );
      _h_Kpip4    = bookHisto1D("h_Kpip4"   ,200,0.,3.2 );
      _h_pipi4    = bookHisto1D("h_pipi4"   ,200,0.,2. );
      _h_Kpi04    = bookHisto1D("h_Kpi04"   ,200,0.,3.2);
      _h_kppim5    = bookHisto1D("h_kppim5"   ,200,0.,3. );
      _h_kppip5    = bookHisto1D("h_kppip5"   ,200,0.,3.1 );
      _h_pippim5   = bookHisto1D("h_pippim5"  ,200,0.,2. );
      _h_kppim6    = bookHisto1D("h_kppim6"   ,200,0.,3.5);
      _h_kppip6    = bookHisto1D("h_kppip6"   ,200,0.,3.5);
      _h_pippim6   = bookHisto1D("h_pippim6"  ,200,0.,2.5);
      _dalitz1     = bookHisto2D("dalitz1"    ,50,0.3,3.2,50,0.3,3.2);
      _dalitz2     = bookHisto2D("dalitz2"    ,50,0.3,3. ,50,0.3,3. );
      _dalitz3     = bookHisto2D("dalitz3"    ,50,0.3,2. ,50,0.07,2. );
      _dalitz4     = bookHisto2D("dalitz4"    ,50,0.3,3.1 ,50,0.07,2. );
      _dalitz5     = bookHisto2D("dalitz5"    ,50,0.,3. ,50,0.,2. );
      _dalitz6     = bookHisto2D("dalitz6"    ,50,0.3,3.5,50,0.07,2.5);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0  ,
			   Particles & Kp  , Particles & Km  , Particles & K0) {
      for(const Particle & p : mother.children()) {
        int id = p.pdgId();
        if ( id == PID::KPLUS ) {
       	  Kp.push_back(p);
	  ++nstable;
	}
	else if (id == PID::KMINUS ) {
	  Km.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIPLUS) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S||id == PID::K0L) {
	  K0.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, Kp , Km, K0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double weight = event.weight();
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::abspid== 411 ||Cuts::abspid== 421 ||Cuts::abspid== 431 )) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, Kp , Km, K0;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km, K0);
	if(meson.pdgId()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	}
	if(abs(meson.pdgId())==421) {
	  if(pim.size()==1&&pip.size()==1&&K0.size()==1) {
	    double mminus = (pim[0].momentum()+K0[0].momentum() ).mass2();
	    double mplus  = (pip[0].momentum()+K0[0].momentum() ).mass2();
	    double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	    _h_plus1  ->fill(mplus,weight);
	    _h_minus1 ->fill(mminus,weight);
	    _h_pipi1  ->fill(mpipi,weight);
	    _dalitz1   ->fill(mplus,mminus,weight);
	  }
	  else if (pip.size()==1&&Km.size()==1&&pi0.size()==1) {
	    double mneut  = (Km[0].momentum()+pip[0].momentum()).mass2();
	    double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pi0[0].momentum()).mass2();
	    _h_neutral2  ->fill(mneut,weight);
	    _h_minus2 ->fill(mminus,weight);
	    _h_pipi2  ->fill(mpipi,weight);
	    _dalitz2->fill(mminus,mneut,weight);
	  }
	}
	else if(abs(meson.pdgId())==411) {
	  if(pip.size()==2&&Km.size()==1) {
	    double mplus  = (Km[0].momentum() +pip[0].momentum()).mass2();
	    double mminus = (Km[0].momentum() +pip[1].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pip[1].momentum()).mass2();
	    if(mplus<mminus) swap(mplus,mminus);
	    _h_Kpilow3 ->fill( mminus,weight);
	    _h_Kpihigh3->fill( mplus,weight);
	    _h_Kpiall3 ->fill( mminus,weight);
	    _h_Kpiall3 ->fill( mplus,weight);
	    _h_pipi3   ->fill( mpipi,weight);
	    _dalitz3->fill(mminus,mpipi, weight);
	  }
	  else if (pip.size()==1&&pi0.size()==1&&K0.size()==1) {
	    double mminus = (K0[0].momentum()+pip[0].momentum()).mass2();
	    double mplus  = (K0[0].momentum()+pi0[0].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pi0[0].momentum()).mass2();
	    _h_Kpip4 ->fill( mminus, weight);
	    _h_pipi4 ->fill( mpipi , weight);
	    _h_Kpi04 ->fill( mplus , weight);
	    _dalitz4->fill(mplus,mpipi, weight);
	  }
	  else if (pim.size()==1&&Kp.size()==1&&pip.size()==1) {
	    double mplus  = (Kp [0].momentum()+pip[0].momentum()).mass2();
	    double mminus = (Kp [0].momentum()+pim[0].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	    _h_kppim5 ->fill(mminus,weight);
	    _h_kppip5 ->fill(mplus ,weight);
	    _h_pippim5->fill(mpipi ,weight);
	    _dalitz5->fill(mminus,mpipi, weight);
	  }
	}
	else if(abs(meson.pdgId())==431) {
	  if (pim.size()==1&&Kp.size()==1&&pip.size()==1) {
	    double mplus  = (Kp [0].momentum()+pip[0].momentum()).mass2();
	    double mminus = (Kp [0].momentum()+pim[0].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	    _h_kppim6 ->fill(mminus,weight);
	    _h_kppip6 ->fill(mplus ,weight);
	    _h_pippim6->fill(mpipi ,weight);
	    _dalitz6->fill(mminus,mpipi, weight);
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_plus1);
      normalize(_h_minus1);
      normalize(_h_pipi1);
      normalize(_dalitz1);
      normalize(_h_minus2);
      normalize(_h_pipi2);
      normalize(_h_neutral2);
      normalize(_dalitz2);
      normalize(_h_Kpilow3);
      normalize(_h_Kpihigh3);
      normalize(_h_Kpiall3);
      normalize(_h_pipi3);
      normalize(_dalitz3);
      normalize(_h_Kpip4);
      normalize(_h_pipi4);
      normalize(_h_Kpi04);
      normalize(_dalitz4);
      normalize(_h_kppim5);
      normalize(_h_kppip5);
      normalize(_h_pippim5);
      normalize(_dalitz5);
      normalize(_h_kppim6);
      normalize(_h_kppip6);
      normalize(_h_pippim6);
      normalize(_dalitz6);
    }
    //@}

    /// @name Histograms
    //@{
    // Histograms for D^0\to \bar{K}^0\pi^+\pi^-
    //m^2_+
    Histo1DPtr _h_plus1;
    //m^2_+
    Histo1DPtr _h_minus1;
    //m^2_{\pi\pi}
    Histo1DPtr _h_pipi1;
    // Dalitz plot
    Histo2DPtr _dalitz1;
    
    // Histograms for D^0\to K^-\pi^+\pi^0
    // Histogram for the K^-\pi^+ mass
    Histo1DPtr _h_minus2;
    // Histogram for the \pi^+\pi^0 mass
    Histo1DPtr _h_pipi2;
    // Histogram for the K^-\pi^0 mass
    Histo1DPtr _h_neutral2;
    // Dalitz plot
    Histo2DPtr _dalitz2;

    // Histograms for D^+\to K^-\pi^+\pi^+
    // Histogram for K^-\pi^+ low
    Histo1DPtr _h_Kpilow3;
    // Histogram for K^-\pi^+ high
    Histo1DPtr _h_Kpihigh3;
    // Histogram for K^-\pi^+ all
    Histo1DPtr _h_Kpiall3;
    // Histogram for \pi^+\pi^-
    Histo1DPtr _h_pipi3;
    // Dalitz plot
    Histo2DPtr _dalitz3;

    // Histograms for D^+\to\bar{K}^0\pi^+\pi^0
    // Histogram for the \bar{K}^0\pi^+ mass
    Histo1DPtr _h_Kpip4;
    // Histogram for the \pi^+\pi^0 mass
    Histo1DPtr _h_pipi4;
    // Histogram for the \bar{K}^0\pi^0 mass
    Histo1DPtr _h_Kpi04;
    // Dalitz plot
    Histo2DPtr _dalitz4;

    // Histograms for D^+\to K^+\pi^-\pi^+
    // Histogram for K^+\pi^-
    Histo1DPtr _h_kppim5;
    // Histogram for K^+\pi^+
    Histo1DPtr _h_kppip5;
    // Histogram for \pi^+\pi^-
    Histo1DPtr _h_pippim5;
    // Dalitz plot
    Histo2DPtr _dalitz5;

    // Histograms for D_s^+\to K^+\pi^-\pi^+
    // Histogram for K^+\pi^-
    Histo1DPtr _h_kppim6;
    // Histogram for K^+\pi^+
    Histo1DPtr _h_kppip6;
    // Histogram for \pi^+\pi^-
    Histo1DPtr _h_pippim6;
    // Dalitz plot
    Histo2DPtr _dalitz6;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_D_Dalitz);


}
