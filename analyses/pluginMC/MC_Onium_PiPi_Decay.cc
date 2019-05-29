// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_Onium_PiPi_Decay : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_Onium_PiPi_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(),"UFS");

    }

    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& pip, Particles& pim,
			   Particles& pi0, Particles & onium) {
      for(const Particle & p : mother.children()) {
        int id = p.pdgId();
      	if ( id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
       	else if (id == PID::PIPLUS) {
       	  pim.push_back(p);
       	  ++nstable;
       	}
       	else if (id == PID::PI0) {
       	  pi0.push_back(p);
       	  ++nstable;
       	}
	else if (abs(id)%1000==443 || abs(id)%1000==553) {
	  onium.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,pip,pim,pi0,onium);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double weight = event.weight();
      // loop over unstable particles
      for(const Particle& vMeson : apply<UnstableParticles>(event, "UFS").particles()) {
	int id = vMeson.pdgId();
	if(id%1000!=443 && id%1000!=553) continue;
	unsigned int nstable(0);
	Particles pip, pim, pi0, onium;
	findDecayProducts(vMeson,nstable,pip,pim,pi0,onium);
	// check for onium
	if(onium.size() !=1 || nstable !=3) continue;
	// check for pipi
	if( ! ((pip.size()==1 && pim.size() ==1) || pi0.size()==2)) continue;
	// check if histos already made
	unsigned int iloc=0; bool found(false);
	while(!found&&iloc<_incoming.size()) {
	  if(_incoming[iloc]==vMeson.pdgId()&&_outgoing[iloc]==onium[0].pdgId()) found=true;
	  else ++iloc;
	}
	// if histos not made, make them
	if(!found) {
	  double twompi = 0.378;
	  double upp    = vMeson.mass()-onium[0].mass();
	  iloc=_incoming.size();
	  _incoming.push_back(vMeson.pdgId());
	  _outgoing.push_back(onium[0].pdgId());
	  ostringstream title;
	  title << "h_" << vMeson.pdgId() << "_" << onium[0].pdgId() << "_";
	  _mpipi.push_back(make_pair(bookHisto1D(title.str()+"mpippim",200,twompi/GeV,upp/GeV),
				     bookHisto1D(title.str()+"mpi0pi0",200,twompi/GeV,upp/GeV)));
	  _hel  .push_back(make_pair(bookHisto1D(title.str()+"hpippim",200,-1.,1.),
				     bookHisto1D(title.str()+"hpi0pi0",200, 0.,1.)));
	}
	// boost to rest frame of the pion pair
	FourMomentum q = vMeson.momentum()-onium[0].momentum();
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(q.betaVec());
	FourMomentum qp = onium[0].momentum();
	FourMomentum ppi= pip.size()==1 ? pip[0].momentum() : pi0[0].momentum();
	qp  = boost.transform(qp);
	ppi = boost.transform(ppi);
	double cX=-ppi.p3().unit().dot(qp.p3().unit());
	if(pi0.size()==2) {
	  _mpipi[iloc].second->fill(q.mass(),weight);
	  _hel  [iloc].second->fill(abs(cX),weight);
	}
	else {
	  _mpipi[iloc].first->fill(q.mass(),weight);
	  _hel  [iloc].first->fill(cX,weight);
	}
      }
    }
    
    
    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize to unity
      for(unsigned int ix=0;ix<_mpipi.size();++ix) {
	normalize(_mpipi[ix].first );
	normalize(_mpipi[ix].second);
	normalize(_hel[ix].first );
	normalize(_hel[ix].second);
      }
    }

    //@}

    /**
     *  Incoming onium states
     */
    vector<long> _incoming;
    
    /**
     *  Outgoing onium states
     */
    vector<long> _outgoing;
    
    /**
     *  Histograms for the \f$\pi^+\pi^-\f$ masses
     */
    vector<pair<Histo1DPtr,Histo1DPtr> > _mpipi;
    
    /**
     *  Histmgrams for the helicity angles
     */
    vector<pair<Histo1DPtr,Histo1DPtr> > _hel;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Onium_PiPi_Decay);


}
