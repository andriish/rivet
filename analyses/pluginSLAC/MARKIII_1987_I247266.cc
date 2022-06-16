// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D -> K pi pi dalitz decays
  class MARKIII_1987_I247266 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKIII_1987_I247266);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      
      book(_h_1_Kmpip,1,1,1);
      book(_h_1_pipi ,1,1,2);
      book(_h_1_Kmpi0,1,1,3);
      book(_dalitz1, "dalitz1",50,0.3,3.2,50,0.,2.);
      
      book(_h_2_K0pip,2,1,1);
      book(_h_2_pipi ,2,1,2);
      book(_h_2_K0pim,2,1,3);
      book(_dalitz2, "dalitz2",50,0.3,3.2,50,0.,2.);
      
      book(_h_3_K0pip,3,1,1);
      book(_h_3_pipi ,3,1,2);
      book(_h_3_K0pi0,3,1,3);
      book(_dalitz3, "dalitz3",50,0.3,3.2,50,0.,2.);
      
      book(_h_4_Kmpip[0],4,1,1);
      book(_h_4_pipi    ,4,1,2);
      book(_h_4_Kmpip[1],4,1,3);
      book(_dalitz4, "dalitz4",50,0.3,2.,50,0.,2.);
      
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0  ,
			   Particles & Kp  , Particles & Km  , Particles & K0) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
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
	else if (id == PID::K0S) {
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
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::abspid== 411 ||Cuts::abspid== 421)) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, Kp , Km, K0;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km, K0);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	}
	if(abs(meson.pid())==421) {
	  if (pip.size()==1&&Km.size()==1&&pi0.size()==1) {
	    double mneut  = (Km[0].momentum()+pip[0].momentum()).mass2();
	    double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pi0[0].momentum()).mass2();
	    _h_1_Kmpip->fill(mneut );
	    _h_1_pipi ->fill(mpipi );
	    _h_1_Kmpi0->fill(mminus);
	    _dalitz1  ->fill(mminus,mpipi);
	  }
	  else if(pim.size()==1&&pip.size()==1&&K0.size()==1) {
	    double mminus = (pim[0].momentum()+K0[0].momentum() ).mass2();
	    double mplus  = (pip[0].momentum()+K0[0].momentum() ).mass2();
	    double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	    _h_2_K0pip->fill(mplus);
	    _h_2_K0pim->fill(mminus);
	    _h_2_pipi ->fill(mpipi);
	    _dalitz2  ->fill(mminus,mpipi); 
	  }
	}
	else if(abs(meson.pid())==411) {
	  if (pip.size()==1&&pi0.size()==1&&K0.size()==1) {
	    double mminus = (K0[0].momentum()+pip[0].momentum()).mass2();
	    double mplus  = (K0[0].momentum()+pi0[0].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pi0[0].momentum()).mass2();
	    _h_3_K0pip->fill( mminus);
	    _h_3_pipi ->fill( mpipi );
	    _h_3_K0pi0->fill( mplus );
	    _dalitz3  ->fill(mplus,mpipi);
	  }
	  else if(pip.size()==2&&Km.size()==1) {
	    double mplus  = (Km[0].momentum() +pip[0].momentum()).mass2();
	    double mminus = (Km[0].momentum() +pip[1].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pip[1].momentum()).mass2();
	    if(mplus<mminus) swap(mplus,mminus);
	    _h_4_Kmpip[1]->fill(mminus);
	    _h_4_Kmpip[0]->fill(mplus );
	    _h_4_pipi    ->fill( mpipi);
	    _dalitz4     ->fill(mminus,mpipi);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_1_Kmpip);
      normalize(_h_1_pipi );
      normalize(_h_1_Kmpi0);
      normalize(_dalitz1  );
      
      normalize(_h_2_K0pip);
      normalize(_h_2_pipi );
      normalize(_h_2_K0pim);
      normalize(_dalitz2  );
      
      normalize(_h_3_K0pip);
      normalize(_h_3_pipi );
      normalize(_h_3_K0pi0);
      normalize(_dalitz3  );
      
      normalize(_h_4_Kmpip[0]);
      normalize(_h_4_pipi    );
      normalize(_h_4_Kmpip[1]);
      normalize(_dalitz4     );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_1_Kmpip, _h_1_pipi, _h_1_Kmpi0;
    Histo2DPtr _dalitz1;
    Histo1DPtr _h_2_K0pip, _h_2_pipi, _h_2_K0pim;
    Histo2DPtr _dalitz2;
    Histo1DPtr _h_3_K0pip, _h_3_pipi, _h_3_K0pi0;
    Histo2DPtr _dalitz3;
    Histo1DPtr _h_4_Kmpip[2], _h_4_pipi;
    Histo2DPtr _dalitz4;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKIII_1987_I247266);

}
