// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_Semi_Leptonic_Decay : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_Semi_Leptonic_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

    }

    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& lp, Particles& lm,
			   Particles& nu, Particles& nub,
			   Particles& out) {
      for(const Particle & p : mother.children()) {
	int id = p.pdgId();
     	if ( id == PID::EMINUS || id == PID::MUON ) {
      	  lm .push_back(p);
      	  ++nstable;
      	}
      	else if (id == PID::EPLUS || id == PID::ANTIMUON) {
      	  lp .push_back(p);
      	  ++nstable;
      	}
	else if ( id == PID::NU_E || id == PID::NU_EBAR ) {
      	  nu .push_back(p);
      	  ++nstable;
      	}
      	else if (id == PID::NU_MU || id == PID::NU_MUBAR ) {
      	  nub.push_back(p);
      	  ++nstable;
      	}
	else if (PID::isMeson(id)) {
	  out.push_back(p);
	  ++nstable;
	}
      	else if ( !p.children().empty() ) {
      	  findDecayProducts(p,nstable,lp,lm,nu,nub,out);
      	}
      	else
      	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double weight = event.weight();
      // loop over unstable particles
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles()) {
	int id = meson.pdgId();
	// spin 0 mesons
	if(!PID::isMeson(id)) continue;
	if(abs(id)%10!=1) continue;
	unsigned int nstable(0);
	Particles lp, lm, nu, nub, out;
	findDecayProducts(meson,nstable,lp,lm,nu,nub,out);
	if(nstable!=3 || out.size()!=1) continue;
	int ilep=0;
	FourMomentum plep,pmnu=out[0].momentum();
	double me2(0.);
	if( lp.size()==1 && nu.size()==1 && out.size()==1 ) {
	  if(nu[0].pdgId()  != -lp[0].pdgId()+1) continue;
	  ilep =  lp[0].pdgId();
	  plep = nu[0].momentum()+lp[0].momentum();
	  pmnu += nu[0].momentum();
	  me2 = lp[0].mass2();
	}
	else if( lm.size()==1 && nub.size()==1 && out.size()==1 ) {
	  if(nub[0].pdgId() != -lm[0].pdgId()-1) continue;
	  ilep =  lm[0].pdgId();
	  plep = nub[0].momentum()+lm[0].momentum();
	  pmnu += nub[0].momentum();
	  me2 = lm[0].mass2();
	}
	else
	  continue;
	// check if histos already exist
	unsigned int iloc=0; bool found(false);
	while(!found&&iloc<_incoming.size()) {
	  if(_incoming[iloc] == id  &&
	     _outgoing[iloc] == out[0].pdgId() &&
	     ilep==_outgoingL[iloc]) found=true; 
	  else ++iloc;
	}
	if(!found) {
	  iloc=_incoming.size();
	  _incoming.push_back(id);
	  _outgoing.push_back(out[0].pdgId());
	  _outgoingL.push_back(ilep);
	  ostringstream title;
	  title << "h_" << abs(id);
	  if(id>0) title << "p";
	  else     title << "m";
	  title << "_" << abs(out[0].pdgId());
	  if(out[0].pdgId()>0) title << "p";
	  else                 title << "m";
	  title << "_" << abs(ilep);
	  if(ilep>0) title << "p";
	  else       title << "m";
	  title << "_";
	  _energy.push_back(bookHisto1D(title.str()+"energy",
					200,0.0,0.5*meson.mass()/MeV));
	  _scale .push_back(bookHisto1D(title.str()+"scale",
					200,0.0,meson.mass()/MeV));
	}
	// add the results to the histogram
	_scale[iloc]->fill(plep.mass()/MeV,weight);
	double ee = 0.5/meson.mass()*(meson.mass2()-pmnu.mass2()+me2);
	_energy[iloc]->fill(ee/MeV,weight);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<_energy.size();++ix) {
	normalize(_energy);
	normalize(_scale );
      }
    }

    //@}


    /// @name Histograms
    //@{
    /**
     *  PDG codes of the decaying mesons
     */ 
    vector<long> _incoming;
    
    /**
     *  PDG codes of the decay products
     */
    vector<long> _outgoing;  
    
    /**
     *  Identidies of the leptons
     */
    vector<long> _outgoingL;
    
    /**
     *  Histograms
     */
    //@{
    /**
     *  The lepton energy
     */
    vector<Histo1DPtr> _energy;
    
    /**
     *  The \f$q\f$ value
     */
    vector<Histo1DPtr> _scale;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Semi_Leptonic_Decay);


}
