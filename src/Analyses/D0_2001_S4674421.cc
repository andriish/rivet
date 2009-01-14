// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Analyses/D0_2001_S4674421.hh" 
#include "AIDA/IDataPoint.h"

namespace Rivet {



  void D0_2001_S4674421::init() {
    _eventsFilledW = 0.0;
    _eventsFilledZ = 0.0;
    _h_dsigdpt_w = bookHistogram1D(1, 1, 1, "$\\mathrm{d}\\sigma / \\mathrm{d}p_\\perp(W)$");
    _h_dsigdpt_z = bookHistogram1D(1, 1, 2, "$\\mathrm{d}\\sigma / \\mathrm{d}p_\\perp(Z)$");

    vector<double> bins(23);
    bins += 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 160, 200;
    _h_dsigdpt_scaled_z = bookHistogram1D("d01-x01-y03", "$\\mathrm{d}\\sigma / \\mathrm{d}(p_\\perp(Z) \\cdot M_W/M_Z)$", bins);
  }



  void D0_2001_S4674421::analyze(const Event& event) {
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
        if (mass/GeV > _mZmin && mass/GeV < _mZmax) {
          Zcount += 1;
          _eventsFilledZ += weight;
          getLog() << Log::DEBUG << "Z #" << Zcount << " pmom.pT() = " << pmom.pT()/GeV << " GeV" << endl;
          _h_dsigdpt_z->fill(pmom.pT()/GeV, weight);
          _h_dsigdpt_scaled_z->fill(pmom.pT()/GeV * _mwmz, weight);
        }
      } else { 
        // There is no Z->ee candidate...
	
        const LeadingParticlesFinalState& enuFS = applyProjection<LeadingParticlesFinalState>(event, "enuFS");
        const LeadingParticlesFinalState& enubFS = applyProjection<LeadingParticlesFinalState>(event, "enubFS"); 
        static size_t Wcount = 0;
	
        // Fill W pT distributions
        /// @todo MERGE THESE TWO BRANCHES BETTER! (including a do-nothing else branch)
        if (enuFS.particles().size() == 2 && enubFS.isEmpty()) {
          const ParticleVector& Wenu = enuFS.particles();
          ParticleVector::const_iterator p = Wenu.begin();
	  
          ParticleVector Wel;
          if (abs(p->pdgId()) == ELECTRON) {
            Wel.push_back(*p);
          }
	  
          FourMomentum pmom = p->momentum();
          ++p;
          pmom += p->momentum();
	  
          if (abs(p->pdgId()) == ELECTRON) {
            Wel.push_back(*p);
          }
    
          ++Wcount;
          _h_dsigdpt_w->fill(pmom.pT()/GeV, weight);
          _eventsFilledW += weight;
        }
        else if (enuFS.isEmpty() && enubFS.particles().size() == 2) {
          const ParticleVector& Wenub = enubFS.particles();
          ParticleVector::const_iterator p = Wenub.begin();
	  
          ParticleVector Wel;
          if (abs(p->pdgId()) == ELECTRON) {
            Wel.push_back(*p); 
          }
	  
          FourMomentum pmom = p->momentum();
          ++p;
          pmom += p->momentum();
	  
          if (abs(p->pdgId()) == ELECTRON) {
            Wel.push_back(*p);
          }
	  
          _h_dsigdpt_w->fill(pmom.pT()/GeV, weight);
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

    // Get W and Z pT integrals
    const double wpt_integral = integral(_h_dsigdpt_w);
    const double zpt_scaled_integral = integral(_h_dsigdpt_scaled_z);

    // Divide and scale ratio histos
    AIDA::IDataPointSet* div = histogramFactory().divide(histoDir() + "/d02-x01-y01", *_h_dsigdpt_w, *_h_dsigdpt_scaled_z); 
    div->setTitle("$[\\mathrm{d}\\sigma/\\mathrm{d}p_\\perp(W)] / [\\mathrm{d}\\sigma/\\mathrm{d}(p_\\perp(Z) \\cdot M_W/M_Z)]$");
    if (xSecW == 0 || wpt_integral == 0 || xSecZ == 0 || zpt_scaled_integral == 0) {
      getLog() << Log::WARN << "Not filling ratio plot because input histos are empty" << endl;
    } else {
      // Scale factor converts event counts to cross-sections, and inverts the
      // branching ratios since only one decay channel has been analysed for each boson.
      const double scalefactor = (xSecW / wpt_integral) / (xSecZ / zpt_scaled_integral) * (_brzee / _brwenu);
      for (int pt = 0; pt < div->size(); ++pt) {
        assert(div->point(pt)->dimension() == 2);
        AIDA::IMeasurement* m = div->point(pt)->coordinate(1);
        m->setValue(m->value() * scalefactor);
        m->setErrorPlus(m->errorPlus() * scalefactor);
        m->setErrorMinus(m->errorPlus() * scalefactor);
      }
    }

    // Normalize non-ratio histos
    normalize(_h_dsigdpt_w, xSecW);
    normalize(_h_dsigdpt_z, xSecZ);
    normalize(_h_dsigdpt_scaled_z, xSecZ);

  }


}
