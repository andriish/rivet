// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_Eta_Decay : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_Eta_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      
      // Book histograms
      double meta[2]={547.45, 957.78};
      for(unsigned int ix=0;ix<2;++ix) {
	std::ostringstream title; title << "_" << ix;
	_mgammagamma .push_back(bookHisto1D("mgammagamma" +title.str(),200,0.,meta[ix]) );
	_mpi0gamma   .push_back(bookHisto1D("mpi0gamma"   +title.str(),200,0.,meta[ix])) ;
	_mpipgamma   .push_back(bookHisto1D("mpipgamma"   +title.str(),200,0.,meta[ix]));
	_mpimgamma   .push_back(bookHisto1D("mpimgamma"   +title.str(),200,0.,meta[ix]));
	_photonenergy.push_back(bookHisto1D("photonenergy"+title.str(),200,0.,meta[ix]));
	_mpippim     .push_back(bookHisto1D("mpippim"     +title.str(),200,0.,meta[ix]));
	_dpippim     .push_back(bookHisto1D("dpippim"     +title.str(),200,200.,meta[ix]));
	_dpi0pi0     .push_back(bookHisto1D("dpi0pi0"     +title.str(),200,200.,meta[ix]));
	_dpi0pip     .push_back(bookHisto1D("dpi0pip"     +title.str(),200,200.,meta[ix]));
	_dpi0pim     .push_back(bookHisto1D("dpi0pim"     +title.str(),200,200.,meta[ix]));
      }
      _dpi0pi0.push_back(bookHisto1D("dpi0pi0_2",200,200.,500.   ));
      _dpippim.push_back(bookHisto1D("dpippim_2",200,200.,500.   ));
      _dpipeta=bookHisto1D("dpipeta",200,500.,meta[1]) ;
      _dpimeta=bookHisto1D("dpimeta",200,500.,meta[1]) ;
      _dpi0eta=bookHisto1D("dpi0eta",200,500.,meta[1]) ;
    }

    void findDecayProducts(const Particle & mother,
                           unsigned int & nstable,
                           Particles& pip, Particles& pim,
			   Particles& pi0, Particles& eta,
			   Particles& gamma) {
      for(const Particle & p : mother.children()) {
	int id = p.pdgId();
        if ( id == PID::ETA ) {
	  eta.push_back(p);
	  ++nstable;
	}
        else if ( id == PID::PHOTON ) {
	  gamma.push_back(p);
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
	else if (id == PID::PI0 ) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S   || id == PID::K0L ||
		 id == PID::KPLUS || id == PID::KMINUS)
	  ++nstable;
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,pip,pim,pi0,eta,gamma);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double weight = event.weight();
      // Loop over f_1 mesons
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::pid==221||Cuts::pid==331)) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, eta, gamma;
	findDecayProducts(meson,nstable,pip, pim, pi0, eta, gamma);
	unsigned int imeson = meson.pdgId()==221 ? 0 : 1;
	// pi0 gamma gamma
	if(nstable==3 && pi0.size()==1 && gamma.size()==2) {
	  _mgammagamma[imeson]->fill((gamma[0].momentum()+gamma[1].momentum()).mass()/MeV,weight);
	  _mpi0gamma  [imeson]->fill((  pi0[0].momentum()+gamma[0].momentum()).mass()/MeV,weight);
	  _mpi0gamma  [imeson]->fill((  pi0[0].momentum()+gamma[1].momentum()).mass()/MeV,weight);
	}  // pi+pi-gamma analysis
	else if(nstable==3 && pip.size()==1 && pim.size()==1 && gamma.size()==1) {
	  FourMomentum ptemp = pip[0].momentum()+pim[0].momentum();
	  double mpipi = ptemp.mass();
	  _mpippim[imeson]->fill(mpipi/MeV,weight);
	  double egamma = 0.5*(meson.mass()*meson.mass()-mpipi*mpipi)/meson.mass();
	  _photonenergy[imeson]->fill(egamma/MeV,weight);
	  _mpipgamma[imeson]->fill((pip[0].momentum()+gamma[0].momentum()).mass()/MeV,weight);
	  _mpimgamma[imeson]->fill((pim[0].momentum()+gamma[0].momentum()).mass()/MeV,weight);
	}
	else if(nstable==3&& pi0.size()==3) {
	  _dpi0pi0[imeson]->fill((pi0[0].momentum()+pi0[1].momentum()).mass()/MeV,weight);
	  _dpi0pi0[imeson]->fill((pi0[0].momentum()+pi0[2].momentum()).mass()/MeV,weight);
	  _dpi0pi0[imeson]->fill((pi0[1].momentum()+pi0[2].momentum()).mass()/MeV,weight);
	}
	else if(nstable==3&& pip.size()==1&&pim.size()==1&&pi0.size()==1) {
	  _dpi0pip[imeson]->fill((pi0[0].momentum()+pip[0].momentum()).mass()/MeV,weight);
	  _dpi0pim[imeson]->fill((pi0[0].momentum()+pim[0].momentum()).mass()/MeV,weight);
	  _dpippim[imeson]->fill((pip[0].momentum()+pim[0].momentum()).mass()/MeV,weight);
	}
	else if(nstable==3&& pi0.size()==2&&eta.size()==1) {
	  _dpi0pi0[2]->fill((pi0[0].momentum()+pi0[1].momentum()).mass()/MeV,weight);
	  _dpi0eta   ->fill((pi0[0].momentum()+eta[0].momentum()).mass()/MeV,weight);
	  _dpi0eta   ->fill((pi0[1].momentum()+eta[0].momentum()).mass()/MeV,weight);
	}
	else if(nstable==3&& pip.size()==1&&pim.size()==1&&eta.size()==1) {
	  _dpippim[2]->fill((pip[0].momentum()+pim[0].momentum()).mass()/MeV,weight);
	  _dpipeta   ->fill((pip[0].momentum()+eta[0].momentum()).mass()/MeV,weight);
	  _dpimeta   ->fill((pim[0].momentum()+eta[0].momentum()).mass()/MeV,weight);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize to unity
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_mgammagamma[ix]);
	normalize(_mpi0gamma[ix]);
	normalize(_mpipgamma[ix]);
	normalize(_mpimgamma[ix]);
	normalize(_mpippim[ix]);
	normalize(_photonenergy[ix]);
	normalize(_dpippim[ix]);
	normalize(_dpi0pi0[ix]);
	normalize(_dpi0pip[ix]);
	normalize(_dpi0pim[ix]);
      }
      normalize(_dpi0pi0[2]);
      normalize(_dpippim[2]);
      normalize(_dpipeta);
      normalize(_dpimeta);
      normalize(_dpi0eta);
    }
    //@}


    /**
     *  Histograms for the decay \f$\eta\to\pi^0\gamma\gamma\f$
     */
    //@{
    /**
     * Histogram for the mass of \f$\gamma\gamma\f$
     */
    vector<Histo1DPtr> _mgammagamma;
    
    /**
     * Histogrma for the mass of \f$\pi^0\gamma\f$
     */
    vector<Histo1DPtr> _mpi0gamma;
    //@}
    
    /**
     *  Histograms for the decay \f$\eta\to\pi^+\pi^-\gamma\f$
     */
    //@{
    /**
     *  Histogram for the mass of \f$\pi^+\gamma\f$
     */
    vector<Histo1DPtr> _mpipgamma;
    
    /**
     *  Histogram for the mass of \f$\pi^-\gamma\f$
     */
    vector<Histo1DPtr> _mpimgamma;
    
    /**
     *  Histogram for the mass of \f$\pi^+\pi^-\f$
     */
    vector<Histo1DPtr> _mpippim;
    
    /**
     *  Histogram for the photon energy
     */
    vector<Histo1DPtr> _photonenergy;
    //@}
    
    /**
     * Histograms for the decay \f$\eta\pi\pi\pi\f$ and \f$\eta'\to\eta\pi\pi\f$.
     */
    //@{
    /**
     *  Histogram for the mass of \f$\pi^+\pi^-\f$
     */
    vector<Histo1DPtr> _dpippim;
    
    /**
     *  Histogram for the mass of \f$\pi^0\pi^0\f$
     */
    vector<Histo1DPtr> _dpi0pi0;
    
    /**
     *  Histogram for the mass of \f$\pi^0\pi^+\f$
     */
    vector<Histo1DPtr> _dpi0pip;
    
    /**
     *  Histogram for the mass of \f$\pi^0\pi^-\f$
     */
    vector<Histo1DPtr> _dpi0pim;
    
    /**
     *  Histogram for the mass of \f$\pi^+\eta\f$
     */
    Histo1DPtr  _dpipeta;
    
    /**
     *  Histogram for the mass of \f$\pi^-\eta\f$
     */
    Histo1DPtr  _dpimeta;
    
    /**
     *  Histogram for the mass of \f$\pi^0\eta\f$
     */
    Histo1DPtr  _dpi0eta;
    //@}
    

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Eta_Decay);


}
