// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Analyses/D0_2001_S4674421.hh" 

namespace Rivet {

  // Booking the histograms:
  void D0_2001_S4674421::init() {	
    Log& log = getLog();
    log << Log::DEBUG << "D0_2001_S4674421::init(): processing" << endl;

    _eventsTried = 0.0;

    _h_dsigdpt_w = bookHistogram1D(1, 1, 1, "dsigma/dpT(W)");
    
    _h_dsigdpt_z = bookHistogram1D(1, 1, 2, "dsigma/dpT(Z)");

    _h_dsigdpt_wz_rat = bookHistogram1D(2, 1, 1, "dsigma/dpT(W) / dsigma/dpT(Z)");


    _eventsFilledW = 0.;
    _eventsFilledZ = 0.;

  }


  void D0_2001_S4674421::analyze(const Event & event) {
      Log& log = getLog();
      log << Log::DEBUG<< "Starting analyzing" << endl;    
 
      
      const double weight = event.weight();
      _eventsTried += weight;
      
      const WZandh& WZbosons = event.applyProjection(_WZproj);


      const ParticleVector& Wens = WZbosons.Wens();

      for (ParticleVector::const_iterator p = Wens.begin(); 
	   p != Wens.end(); ++p) {
        FourMomentum pmom(p->getMomentum());
	_h_dsigdpt_w->fill(pmom.polarRadius(), weight);
	_eventsFilledW += weight;
      }

      
      const ParticleVector& Zees = WZbosons.Zees();

      for (ParticleVector::const_iterator p = Zees.begin(); 
	   p != Zees.end(); ++p) {
        FourMomentum pmom(p->getMomentum());
	static int count =0;
	cout << " Z #" << count++ << "   pmom.polarRadius()=" << pmom.polarRadius() << endl;
	_h_dsigdpt_z->fill(pmom.polarRadius(), weight);
	_eventsFilledZ += weight;
      }


      log << Log::DEBUG << "Finished analyzing" << endl;
}


void D0_2001_S4674421::finalize() { 

  double xSecPerEvent = crossSection() / _eventsTried;
  // HepData data is in nb, crossSection returns pb.
  /// @todo Choose consistent units set
  xSecPerEvent = 0.001 * xSecPerEvent;

  const double xSecW = xSecPerEvent * _eventsFilledW;
  const size_t nBinsW = _h_dsigdpt_w->axis().bins();
  double hAreaW = 0.0;
  for (size_t iBin = 0; iBin != nBinsW; ++iBin) {
    hAreaW += _h_dsigdpt_w->binHeight(iBin) * _h_dsigdpt_w->axis().binWidth(iBin);
  }
  _h_dsigdpt_w->scale(xSecW / hAreaW);


  const double xSecZ = xSecPerEvent * _eventsFilledZ;
  const size_t nBinsZ = _h_dsigdpt_z->axis().bins();
  double hAreaZ = 0.0;
  for (size_t iBin = 0; iBin != nBinsZ; ++iBin) {
    hAreaZ += _h_dsigdpt_z->binHeight(iBin) * _h_dsigdpt_z->axis().binWidth(iBin);
  }  
  _h_dsigdpt_z->scale(xSecZ / hAreaZ);


  _h_dsigdpt_wz_rat =  histogramFactory().divide("/D0_2001_S4674421/d02-x01-y01", *_h_dsigdpt_w, *_h_dsigdpt_z);

  _h_dsigdpt_wz_rat->scale(_mwmz * _brwenu * _brzee);

    
  

  }

}
