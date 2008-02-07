// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Analyses/D0_2001_S4674421.hh" 

namespace Rivet {


  void D0_2001_S4674421::init() {	
    Log& log = getLog();
    log << Log::TRACE << "D0_2001_S4674421::init(): processing" << endl;
    _eventsFilledW = 0.0;
    _eventsFilledZ = 0.0;
    _h_dsigdpt_w = bookHistogram1D(1, 1, 1, "dsigma/dpT(W)");
    _h_dsigdpt_z = bookHistogram1D(1, 1, 2, "dsigma/dpT(Z)");
    _h_dsigdpt_wz_rat = bookHistogram1D(2, 1, 1, "dsigma/dpT(W) / dsigma/dpT(Z)");
  }


  void D0_2001_S4674421::analyze(const Event & event) {
      Log& log = getLog();
      log << Log::DEBUG << "Starting analyzing" << endl;    

      const double weight = event.weight();
      const WZandh& WZbosons = event.applyProjection(_WZproj);

      // Fill W pT distributions
      const ParticleVector& Wens = WZbosons.Wens();
      for (ParticleVector::const_iterator p = Wens.begin(); p != Wens.end(); ++p) {
        FourMomentum pmom = p->getMomentum();
        _h_dsigdpt_w->fill(pmom.pT(), weight);
        _eventsFilledW += weight;
      }

      // Fill Z pT distributions
      size_t Zcount = 0;      
      const ParticleVector& Zees = WZbosons.Zees();
      for (ParticleVector::const_iterator p = Zees.begin(); p != Zees.end(); ++p) {
        FourMomentum pmom = p->getMomentum();
        log << Log::DEBUG << "Z #" << ++Zcount << " pmom.pT() = " << pmom.pT() << endl;
        _h_dsigdpt_z->fill(pmom.pT(), weight);
        _eventsFilledZ += weight;
      }

      log << Log::DEBUG << "Finished analyzing" << endl;
  }


  void D0_2001_S4674421::finalize() { 
    // Get cross-section per event (i.e. per unit weight) from generator
    const double xSecPerEvent = crossSection()/nanobarn / sumOfWeights();

    // Correct W pT distribution to W cross-section
    const double xSecW = xSecPerEvent * _eventsFilledW;
    cout << "xSecW=" << xSecW << endl;

    // Correct Z pT distribution to Z cross-section
    const double xSecZ = xSecPerEvent * _eventsFilledZ;
    cout << "xSecZ=" << xSecZ << endl;

    _h_dsigdpt_wz_rat = histogramFactory().divide("/D0_2001_S4674421/d02-x01-y01", *_h_dsigdpt_w, *_h_dsigdpt_z);
    _h_dsigdpt_wz_rat->scale(xSecW/xSecZ * _mwmz * _brzee / _brwenu);

    normalize(_h_dsigdpt_w, xSecW);

    normalize(_h_dsigdpt_z, xSecZ);

  }

}
