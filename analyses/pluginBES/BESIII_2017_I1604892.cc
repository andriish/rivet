// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief  e+ e- -> pi+ pi- J/psi
  class BESIII_2017_I1604892 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1604892);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      unsigned int iloc=1;
      if      (isCompatibleWithSqrtS(4.23)) iloc=1;
      else if (isCompatibleWithSqrtS(4.26)) iloc=2;
      for(unsigned int ix=0;ix<2;++ix)
      	book(_h_mass[ix],1,iloc,1+ix);
      book(_h_cTheta,2,1,1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles Jpsi,other;
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
	for (auto & psi : Jpsi) {
	  // J/psi already in list 
	  if (fuzzyEquals(psi.momentum(),parent.momentum())) {
	    found=true;
	    break;
	  }
	}
	if(!found) Jpsi.push_back(parent);
      }
      if(Jpsi.size()!=1 || other.size()!=2) vetoEvent;
      if(other[0].pid()!=-other[1].pid() || other[0].abspid()!=PID::PIPLUS) vetoEvent;
      _h_mass[0]->fill((other[0].momentum()+other[1].momentum()).mass());
      for(unsigned int ix=0;ix<2;++ix) {
	FourMomentum pZ = Jpsi[0].momentum()+other[ix].momentum(); 
	double mZ = pZ.mass();
	_h_mass[1]->fill(mZ);
	if(mZ>3.86 && mZ<3.92) {
	  LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(pZ.betaVec());
	  Vector3 axis1 = pZ.p3().unit();
	  Vector3 axis2 = boost.transform(other[ix].momentum()).p3().unit();
	  _h_cTheta->fill(abs(axis1.dot(axis2)));
	}
      }
  }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h_mass[ix],1.,false);
      normalize(_h_cTheta,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2];
    Histo1DPtr _h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1604892);

}
// BEGIN YODA_SCATTER2D_V2 /REF/BESIII_2017_I1604892/d01-x01-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BESIII_2017_I1604892/d01-x01-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BESIII_2017_I1604892/d01-x02-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BESIII_2017_I1604892/d01-x02-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BESIII_2017_I1604892/d02-x01-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BESIII_2017_I1604892/d02-x01-y02
