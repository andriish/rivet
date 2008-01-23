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
    // Apply cross-section normaisation corrections to distributions
    /// @todo For version 2, we need some more automated/built-in system for applying these corrections

    // Get cross-section per event (i.e. per unit weight) from generator
    /// @todo Make units of xSecPerEvent explicit i.e. use crossSection()/GeV or crossSection()/MeV
    const double xSecPerEvent = crossSection()/1000.0 / sumOfWeights();

    // Correct W pT distribution to W cross-section
    const double xSecW = xSecPerEvent * _eventsFilledW;
    const size_t nBinsW = _h_dsigdpt_w->axis().bins();
    double hAreaW = 0.0;
    for (size_t iBin = 0; iBin != nBinsW; ++iBin) {
      hAreaW += _h_dsigdpt_w->binHeight(iBin) * _h_dsigdpt_w->axis().binWidth(iBin);
    }
    _h_dsigdpt_w->scale(xSecW / hAreaW);

    // Correct Z pT distribution to Z cross-section
    const double xSecZ = xSecPerEvent * _eventsFilledZ;
    const size_t nBinsZ = _h_dsigdpt_z->axis().bins();
    double hAreaZ = 0.0;
    for (size_t iBin = 0; iBin != nBinsZ; ++iBin) {
      hAreaZ += _h_dsigdpt_z->binHeight(iBin) * _h_dsigdpt_z->axis().binWidth(iBin);
    }
    _h_dsigdpt_z->scale(xSecZ / hAreaZ);

    // Make pT_W/pT_Z ratio histogram
    _h_dsigdpt_wz_rat = histogramFactory().divide("/D0_2001_S4674421/d02-x01-y01", *_h_dsigdpt_w, *_h_dsigdpt_z);
    _h_dsigdpt_wz_rat->scale(_mwmz * _brwenu * _brzee);
  }

}
