// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief e+ e- -> e+e- / mu+mu-
  class BESIII_2022_I2102287 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2102287);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_c[ix],"TMP/c_"+toString(ix+1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const double cos34 = cos(34./180.*M_PI);
      // get the axis, direction of incoming positron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis1 = beams.first .momentum().p3().unit();
      Vector3 axis2 = beams.second.momentum().p3().unit();
      if(beams.first.pid()<0) swap(axis1,axis2);
      _c[2]->fill((beams.first.momentum()+beams.second.momentum()).mass());
      // loop over FS particles
      Particles fs = apply<FinalState>(event,"FS").particles();
      Particles em,ep,mm,mp;
      for(const Particle & p : fs) {
	if(p.abspid()==PID::MUON) {
	  if(p.pid()>0) mm.push_back(p);
	  else          mp.push_back(p);
	}
	else if(p.abspid()==PID::ELECTRON) {
	  if(abs(axis1.dot(p.p3().unit()))>cos34) continue;
	  if(p.pid()>0) em.push_back(p);
	  else          ep.push_back(p);
	}
	else if(p.pid()!=PID::GAMMA)
	  vetoEvent;
      }
      if(em.size()==1 && ep.size()==1 && mm.size()==0 && mp.size()==0) {
	_c[0]->fill();
      }
      else if(mm.size()==1 && mp.size()==1 && em.size()==0 && ep.size()==0) {
	_c[1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double energy = _c[2]->val()/ sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) {
	double sig_m = _c[ix]  ->val()*fact;
	double err_m = _c[ix]  ->err()*fact;
	Scatter2D temphisto(refData(1+ix, 1, 1));
	Scatter2DPtr cross;
	book(cross, 1+ix, 1, 1);
	double deltaE=1e30;
	unsigned int ipoint=100000000;
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  double test = abs(temphisto.point(b).x()-energy);
	  if(test<deltaE) {
	    deltaE=test;
	    ipoint=b;
	  }
	}
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  if (b!=ipoint)
	    cross  ->addPoint(x, 0., ex, make_pair(0.,.0));
	  else
	    cross  ->addPoint(x, sig_m, ex, make_pair(err_m,err_m));
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2102287);

}
