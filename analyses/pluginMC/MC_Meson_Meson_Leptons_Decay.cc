// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_Meson_Meson_Leptons_Decay : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_Meson_Meson_Leptons_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      
      // Initialise and register projections
      declare(UnstableParticles(),"UFS");

    }

    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& lp, Particles& lm,
			   Particles& scalar,
			   Particles& vector) {
      for(const Particle & p : mother.children()) {
        int id = p.pdgId();
	if ( id == PID::EMINUS || id == PID::MUON ) {
       	  lm.push_back(p);
       	  ++nstable;
       	}
	else if (id == PID::EPLUS || id == PID::ANTIMUON) {
	  lp.push_back(p);
	  ++nstable;
	}
	else if (abs(id)%10==1 && PID::isMeson(id)) {
	  scalar.push_back(p);
	  ++nstable;
	}
	else if ((abs(id)%10==3 && PID::isMeson(id)) ||
		 id==PID::PHOTON ) {
	  vector.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,lp,lm,scalar,vector);
	}
	else
	  ++nstable;
      }
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& iMeson : apply<UnstableParticles>(event, "UFS").particles()) {
	// only consider scalar/vector mesons
	long pid = iMeson.pdgId();
	if(!PID::isMeson(pid)) continue;
	if(abs(pid)%10!=3 and abs(pid)%10!=1 ) continue;
	Particles lp,lm,scalar,vector;
	unsigned int nstable(0);
	findDecayProducts(iMeson,nstable,lp,lm,scalar,vector);
	if(nstable!=3 || lp.size()!=1 || lm.size()!=1 || lp[0].pdgId()!=-lm[0].pdgId()) continue;
	if(scalar.size()==1) {
	  // check if we already have this decay
	  unsigned int ix=0; bool found(false); 
	  while(!found&&ix<_incomingV.size()) {
	    if(_incomingV[ix]==pid && _outgoingP[ix]==scalar[0].pdgId() &&
	       _outgoingf_V[ix]==lm[0].pdgId()) {
	      found=true;
	    }
	    else {
	      ++ix;
	    }
	  }
	  // create a new graph if needed
	  if(!found) {
	    ix=_incomingV.size();
	    _incomingV.push_back(pid);
	    _outgoingP.push_back(scalar[0].pdgId());
	    _outgoingf_V.push_back(lm[0].pdgId());
	    ostringstream title;
	    title << "h_" << abs(pid);
	    if(pid>0) title << "p";
	    else      title << "m";
	    title << "_" << abs(scalar[0].pdgId());
	    if(scalar[0].pdgId()>0) title << "p";
	    else                    title << "m";
	    title << "_" << lm[0].pdgId() << "_";
	    _mff_V   .push_back(bookHisto1D(title.str()+"mff"   , 200, 0., iMeson.mass()));
	    _mPf   .push_back(bookHisto1D(title.str()+"mPf"   , 200, 0., iMeson.mass()));
	    _mPfbar.push_back(bookHisto1D(title.str()+"mPfbar", 200, 0., iMeson.mass()));
	  }
	  // add the results to the histogram
	  _mff_V   [ix]->fill((lm    [0].momentum()+lp[0].momentum()).mass(),event.weight());
	  _mPf   [ix]->fill((scalar[0].momentum()+lm[0].momentum()).mass(),event.weight());
	  _mPfbar[ix]->fill((scalar[0].momentum()+lp[0].momentum()).mass(),event.weight());
	}
	else if(vector.size()==1) {
	  // check if we already have this decay
	  unsigned int ix=0; bool found(false); 
	  while(!found&&ix<_incoming_P.size()) {
	    if(_incoming_P[ix]==pid && _outgoingV[ix]==vector[0].pdgId() &&
	       _outgoingf_P[ix]==lm[0].pdgId()) {
	      found=true;
	    }
	    else {
	      ++ix;
	    }
	  }
	  // create a new graph if needed
	  if(!found) {
	    ix=_incoming_P.size();
	    _incoming_P.push_back(pid);
	    _outgoingV.push_back(vector[0].pdgId());
	    _outgoingf_P.push_back(lm[0].pdgId());
	    ostringstream title;
	    title << "h2_" << abs(pid);
	    if(pid>0) title << "p";
	    else      title << "m";
	    title << "_" << abs(vector[0].pdgId());
	    if(vector[0].pdgId()>0) title << "p";
	    else                    title << "m";
	    title << "_" << lm[0].pdgId() << "_";
	    _mff_P   .push_back(bookHisto1D(title.str()+"mff"   , 200, 0., iMeson.mass()));
	    _mVf   .push_back(bookHisto1D(title.str()+"mVf"   , 200, 0., iMeson.mass()));
	    _mVfbar.push_back(bookHisto1D(title.str()+"mVfbar", 200, 0., iMeson.mass()));
	  }
	  // add the results to the histogram
	  _mff_P   [ix]->fill((lm    [0].momentum()+lp[0].momentum()).mass(),event.weight());
	  _mVf   [ix]->fill((vector[0].momentum()+lm[0].momentum()).mass(),event.weight());
	  _mVfbar[ix]->fill((vector[0].momentum()+lp[0].momentum()).mass(),event.weight());
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      
      // normalize to unity V->P
      for(unsigned int ix=0;ix<_mff_V.size();++ix) {
	normalize(_mff_V);
	normalize(_mPf);
	normalize(_mPfbar);
      }
      // normalize to unity P->V
      for(unsigned int ix=0;ix<_mff_P.size();++ix) {
	normalize(_mff_P);
	normalize(_mVf);
	normalize(_mVfbar);
      }

    }

    //@}



    /// @name Histograms for V -> P
    //@{
    /**
     *  PDG codes of the incoming particles
     */
    vector<long> _incomingV;
    
    /**
     *  PDG codes of the outgoing pseudoscalar mesons
     */
    vector<long> _outgoingP;
    
    /**
     *  PDG codes of the outgoing fermion
     */
    vector<long> _outgoingf_V;
    
    /**
     *  Histograms for the mass of the fermion-antifermion pair
     */
    vector<Histo1DPtr> _mff_V;
    
    /**
     *  Histograms for the masses of the pseudoscalar and the fermion
     */
    vector<Histo1DPtr> _mPf;
    
    /**
     *  Histograms for the masses of the pseudoscalar and the antifermion
     */
    vector<Histo1DPtr> _mPfbar;
    //@}
    
    /// @name Histograms P->V
    //@{
    /**
     *  PDG codes of the incoming_P particles
     */
    vector<long> _incoming_P;
    
    /**
     *  PDG codes of the outgoing vector mesons
     */
    vector<long> _outgoingV;
    
    /**
     *  PDG codes of the outgoing fermion
     */
    vector<long> _outgoingf_P;
    
    /**
     *  Histograms for the mass of the fermion-antifermion pair
     */
    vector<Histo1DPtr> _mff_P;
    
    /**
     *  Histograms for the masses of the vector and the fermion
     */
    vector<Histo1DPtr> _mVf;
    
    /**
     *  Histograms for the masses of the vector and the antifermion
     */
    vector<Histo1DPtr> _mVfbar;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Meson_Meson_Leptons_Decay);


}
