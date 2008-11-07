// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Analyses/D0_2001_S4674421.hh" 
#include "AIDA/IDataPoint.h"

namespace Rivet {


  void D0_2001_S4674421::init() {	
    getLog() << Log::TRACE << "D0_2001_S4674421::init(): processing" << endl;
    _eventsFilledW = 0.0;
    _eventsFilledZ = 0.0;
    _h_dsigdpt_w = bookHistogram1D(1, 1, 1, "$\\mathrm{d}\\sigma / \\mathrm{d}p_\\perp(W)$");
    _h_dsigdpt_z = bookHistogram1D(1, 1, 2, "$\\mathrm{d}\\sigma / \\mathrm{d}p_\\perp(Z)$");

    /// @todo Use autobooking
    const double binning[] = {0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 
                              25., 30., 35., 40., 50., 60., 70., 80., 100., 120., 160., 200.};
    vector<double> bins(23);
    for (size_t i = 0; i < bins.size(); ++i) {
      bins[i] = binning[i];
    }
    _h_dsigdpt_scaled_z = bookHistogram1D("d01-x01-y03","$\\mathrm{d}\\sigma / \\mathrm{d}(p_\\perp(Z) \\cdot M_W/M_Z)$", bins);
    /// @todo Remove need for temp histo title with YODA
    _h_temp = bookHistogram1D("temp", "$[\\mathrm{d}\\sigma/\\mathrm{d}p_\\perp(W)] / [\\mathrm{d}\\sigma/\\mathrm{d}(p_\\perp(Z) \\cdot M_W/M_Z)]$", bins);
  }


  void D0_2001_S4674421::analyze(const Event & event) {
      const double weight = event.weight();
      
      const LeadingParticlesFinalState& eeFS = applyProjection<LeadingParticlesFinalState>(event, "eeFS");
      
      // If there is a Z candidate
      if (eeFS.particles().size() == 2) {
        static size_t Zcount = 0;
        // Fill Z pT distributions
        const ParticleVector& Zees = eeFS.particles();

        ParticleVector::const_iterator p = Zees.begin();
        FourMomentum pmom = p->momentum();
        p++;
        pmom += p->momentum();
        double mass = sqrt(pmom.invariant());
        if (mass > _mZmin && mass < _mZmax) {
          Zcount += 1;
          _eventsFilledZ += weight;
          getLog() << Log::DEBUG << "Z #" << Zcount << " pmom.pT() = " << pmom.pT() << endl;
          _h_dsigdpt_z->fill(pmom.pT(), weight);
          _h_dsigdpt_scaled_z->fill(pmom.pT()*_mwmz, weight);
        }
      } else { 
        // There is no Z->ee candidate...
	
        const LeadingParticlesFinalState& enuFS = applyProjection<LeadingParticlesFinalState>(event, "enuFS");
        const LeadingParticlesFinalState& enubFS = applyProjection<LeadingParticlesFinalState>(event, "enubFS"); 
        static size_t Wcount = 0;
	
        // Fill W pT distributions
        if (enuFS.particles().size() == 2 && enubFS.isEmpty()) {
          const ParticleVector& Wenu = enuFS.particles();
          ParticleVector::const_iterator p = Wenu.begin();
	  
          ParticleVector Wel;
          if (abs(p->pdgId()) == 11) Wel.push_back(*p);
	  
          FourMomentum pmom = p->momentum();
          ++p;
          pmom += p->momentum();
	  
          if (abs(p->pdgId()) == 11) Wel.push_back(*p);
    
          ++Wcount;
          _h_dsigdpt_w->fill(pmom.pT(), weight);
          _eventsFilledW += weight;
        }
        else if (enuFS.isEmpty() && enubFS.particles().size() == 2) {
          const ParticleVector& Wenub = enubFS.particles();
          ParticleVector::const_iterator p = Wenub.begin();
	  
          ParticleVector Wel;
          if (abs(p->pdgId()) == 11) Wel.push_back(*p);
	  
          FourMomentum pmom = p->momentum();
          ++p;
          pmom += p->momentum();
	  
          if (abs(p->pdgId()) == 11) Wel.push_back(*p);
	  
          _h_dsigdpt_w->fill(pmom.pT(), weight);
          _eventsFilledW += weight;
        }
      }

  }


  void D0_2001_S4674421::finalize() { 
    // Get cross-section per event (i.e. per unit weight) from generator
    const double xSecPerEvent = crossSection()/picobarn / sumOfWeights();
    
    // Correct W pT distribution to W cross-section
    const double xSecW = xSecPerEvent * _eventsFilledW;
    
    // Correct Z pT distribution to Z cross-section
    const double xSecZ = xSecPerEvent * _eventsFilledZ;

    /// @todo This should be simpler with YODA!
    _h_temp = histogramFactory().divide(getName() + "/temp", *_h_dsigdpt_w, *_h_dsigdpt_scaled_z);    
    const double wpt_integral = integral(_h_dsigdpt_w);
    normalize(_h_dsigdpt_w, xSecW);
    normalize(_h_dsigdpt_z, xSecZ);
    const double zpt_scaled_integral = integral(_h_dsigdpt_scaled_z);
    normalize(_h_dsigdpt_scaled_z, xSecZ);    
    _h_temp->scale( (xSecW / wpt_integral) / (xSecZ / zpt_scaled_integral)  *  _brzee / _brwenu);
    std::vector<double> _x, _y, _ex, _ey;
    for ( int i = 0, N = _h_temp->axis().bins(); i < N; ++i ) {
      _x.push_back((_h_temp->axis().binLowerEdge(i) + _h_temp->axis().binUpperEdge(i))/2.0);
      _ex.push_back(_h_temp->axis().binWidth(i)/2.0);
      _y.push_back(_h_temp->binHeight(i)); 
      _ey.push_back(_h_temp->binError(i)); 
    }
    _dset_dsigpt_wz_rat = datapointsetFactory().createXY(getName() + "/d02-x01-y01", _h_temp->title(), _x, _y, _ex, _ey);
  }


}
