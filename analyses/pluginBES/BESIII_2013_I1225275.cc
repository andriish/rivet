// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief 
  class BESIII_2013_I1225275 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1225275);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      // histograms
      book(_h_sigma, 1,1,1);
      for(unsigned int ix=0;ix<3;++ix)
	book(_h_mass[ix],2,1,1+ix);
      book(_h_mass[3],3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles PSI,other;
      for(const Particle & p : fs) {
	Particle parent=p;
       	while(!parent.parents().empty()) {
      	  parent=parent.parents()[0];
	  if(parent.abspid()==PID::JPSI) break;
	}
	if(parent.abspid()!=PID::JPSI) {
	  other.push_back(p);
	  continue;
	}
	bool found=false;
	for (auto & p : PSI) {
      	  // psi already in list 
      	  if (fuzzyEquals(p.momentum(),parent.momentum())) {
      	    found=true;
      	    break;
      	  }
      	}
      	if(!found) PSI.push_back(parent);
      }
      // J/psi + 2 other particles
      if(PSI.size()!=1 || other.size()!=2) vetoEvent;
      // other particles pi+ pi-
      if(!(other[0].pid()==-other[1].pid() && other[0].abspid()==PID::PIPLUS)) vetoEvent;
      if(other[0].pid()<0) swap(other[0],other[1]);
      // cross secrtion
      _h_sigma->fill(4.26);
      // fill the mass plots
      double mJpsipi[2] = {(other[0].momentum()+PSI[0].momentum()).mass(),
			   (other[1].momentum()+PSI[0].momentum()).mass()};
      for(unsigned int ix=0;ix<2;++ix) _h_mass[ix]->fill(mJpsipi[ix]);
      _h_mass[2]->fill((other[0].momentum()+other[1].momentum()).mass());
      _h_mass[3]->fill(max(mJpsipi[0],mJpsipi[1]));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_sigma, crossSection()/ sumOfWeights() /picobarn); 
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h_mass[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_sigma,_h_mass[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1225275);

}
