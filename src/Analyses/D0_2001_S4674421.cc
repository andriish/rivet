// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Analyses/D0_2001_S4674421.hh" 
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


    //  - @c _mwmz = ratio of \f$ mW/mZ \f$ used in the publication analysis
    //  - @c _brwenu = ratio of \f$ BR(W->e,nu) \f$ used in the publication analysis
    //  - @c _brzee = ratio of \f$ BR(Z->ee) \f$ used in the publication analysis
    //  - @c _mZmin = lower Z mass cut used in the publication analysis
    //  - @c _mZmax = upper Z mass cut used in the publication analysis
    D0_2001_S4674421::D0_2001_S4674421()
      : Analysis("D0_2001_S4674421"),
        _mwmz(0.8820), _brwenu(0.1073), _brzee(0.033632), 
        _mZmin(75.*GeV), _mZmax(105.*GeV)
    { 

      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
      //const FinalState fs(-3.0, 3.0); 
      FinalState fs(-5.0, 5.0); //corrected for detector acceptance
      addProjection(fs, "FS");

      // Z -> e- e+
      LeadingParticlesFinalState eeFS(fs, -2.5, 2.5, 0.); //20.);
      eeFS.addParticleIdPair(ELECTRON);
      addProjection(eeFS, "eeFS");
      
      // W- -> e- nu_e~
      LeadingParticlesFinalState enuFS(fs, -2.5, 2.5, 0.); //25.);
      enuFS.addParticleId(ELECTRON).addParticleId(NU_EBAR);
      addProjection(enuFS, "enuFS");
      
      // W+ -> e+ nu_e
      LeadingParticlesFinalState enubFS(fs, -2.5, 2.5, 0.); //25.);
      enubFS.addParticleId(POSITRON).addParticleId(NU_E);
      addProjection(enubFS, "enubFS");

      // Remove neutrinos for isolation of final state particles
      VetoedFinalState vfs(fs);
      vfs.vetoNeutrinos();
      addProjection(vfs, "VFS");

    }    




  void D0_2001_S4674421::init() {
    _eventsFilledW = 0.0;
    _eventsFilledZ = 0.0;
    _h_dsigdpt_w = 
      bookHistogram1D(1, 1, 1, "$\\mathrm{d}{\\sigma} / \\mathrm{d}{p_\\perp(W)}$", 
                      "$p_\\perp$ / GeV/$c$", "$\\mathrm{d}{\\sigma}/\\mathrm{d}{p_\\perp(W)}$");
    _h_dsigdpt_z = 
      bookHistogram1D(1, 1, 2, "$\\mathrm{d}{\\sigma} / \\mathrm{d}{p_\\perp(Z)}$", 
                      "$p_\\perp$ / GeV/$c$", "$\\mathrm{d}{\\sigma}/\\mathrm{d}{p_\\perp(Z)}$");

    vector<double> bins(23);
    bins += 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 160, 200;
    _h_dsigdpt_scaled_z = 
      bookHistogram1D("d01-x01-y03", "$\\mathrm{d}{\\sigma} / \\mathrm{d}{(p_\\perp(Z) \\cdot M_W/M_Z)}$", 
                      "$p_\\perp(W)/M_W$", "$R_p_\\perp$", bins);
  }



  void D0_2001_S4674421::analyze(const Event& event) {
    const double weight = event.weight();
    
    const LeadingParticlesFinalState& eeFS = applyProjection<LeadingParticlesFinalState>(event, "eeFS");
    if (eeFS.particles().size() == 2) {
      // If there is a Z candidate:
      static size_t Zcount = 0;
      // Fill Z pT distributions
      const ParticleVector& Zdaughters = eeFS.particles();
      const FourMomentum pmom = Zdaughters[0].momentum() + Zdaughters[1].momentum();
      double mass = sqrt(pmom.invariant());
      if (mass/GeV > _mZmin && mass/GeV < _mZmax) {
        ++Zcount;
        _eventsFilledZ += weight;
        getLog() << Log::DEBUG << "Z #" << Zcount << " pmom.pT() = " << pmom.pT()/GeV << " GeV" << endl;
        _h_dsigdpt_z->fill(pmom.pT()/GeV, weight);
        _h_dsigdpt_scaled_z->fill(pmom.pT()/GeV * _mwmz, weight);
      }
    } else { 
      // There is no Z -> ee candidate... so this must be a W event, right?
      const LeadingParticlesFinalState& enuFS = applyProjection<LeadingParticlesFinalState>(event, "enuFS");
      const LeadingParticlesFinalState& enubFS = applyProjection<LeadingParticlesFinalState>(event, "enubFS"); 
      static size_t Wcount = 0;
      
      // Fill W pT distributions
      ParticleVector Wdaughters;
      if (enuFS.particles().size() == 2 && enubFS.isEmpty()) {
        Wdaughters = enuFS.particles();
      } else if (enuFS.isEmpty() && enubFS.particles().size() == 2) {
        Wdaughters = enubFS.particles();
      }
      if (! Wdaughters.empty()) {
        assert(Wdaughters.size() == 2);
        const FourMomentum pmom = Wdaughters[0].momentum() + Wdaughters[1].momentum();
        ++Wcount;
        _eventsFilledW += weight;
        _h_dsigdpt_w->fill(pmom.pT()/GeV, weight);
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
