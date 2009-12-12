// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief D0 Run I differential W/Z boson cross-section analysis
  /// @author Lars Sonnenschein
  class D0_2001_S4674421 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor.
    D0_2001_S4674421() : Analysis("D0_2001_S4674421") {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
    }
 
 
    /// @name Analysis methods
    //@{
 
    void init() {
      // Final state projection
      FinalState fs(-5.0, 5.0); // corrected for detector acceptance
      addProjection(fs, "FS");

      // Z -> e- e+
      LeadingParticlesFinalState eeFS(FinalState(-2.5, 2.5, 0.)); //20.);
      eeFS.addParticleIdPair(ELECTRON);
      addProjection(eeFS, "eeFS");
   
      // W- -> e- nu_e~
      LeadingParticlesFinalState enuFS(FinalState(-2.5, 2.5, 0.)); //25.);
      enuFS.addParticleId(ELECTRON).addParticleId(NU_EBAR);
      addProjection(enuFS, "enuFS");
   
      // W+ -> e+ nu_e
      LeadingParticlesFinalState enubFS(FinalState(-2.5, 2.5, 0.)); //25.);
      enubFS.addParticleId(POSITRON).addParticleId(NU_E);
      addProjection(enubFS, "enubFS");

      // Remove neutrinos for isolation of final state particles
      VetoedFinalState vfs(fs);
      vfs.vetoNeutrinos();
      addProjection(vfs, "VFS");

      // Counters
      _eventsFilledW = 0.0;
      _eventsFilledZ = 0.0;

      // Histograms
      _h_dsigdpt_w = bookHistogram1D(1, 1, 1);
      _h_dsigdpt_z = bookHistogram1D(1, 1, 2);
      /// @todo Can't take this from ref data? 
      vector<double> bins(23);
      bins += 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 160, 200;
      _h_dsigdpt_scaled_z = bookHistogram1D("d01-x01-y03", bins);
    }



    void analyze(const Event& event) {
      const double weight = event.weight();

      const LeadingParticlesFinalState& eeFS = applyProjection<LeadingParticlesFinalState>(event, "eeFS");
      if (eeFS.particles().size() == 2) {
        // If there is a Z candidate:
        static size_t Zcount = 0;
        // Fill Z pT distributions
        const ParticleVector& Zdaughters = eeFS.particles();
        const FourMomentum pmom = Zdaughters[0].momentum() + Zdaughters[1].momentum();
        double mass = sqrt(pmom.invariant());
        if (inRange(mass/GeV, 75.0, 105.0)) {
          ++Zcount;
          _eventsFilledZ += weight;
          getLog() << Log::DEBUG << "Z #" << Zcount << " pmom.pT() = " << pmom.pT()/GeV << " GeV" << endl;
          const double MW_MZ = 0.8820; // Ratio M_W/M_Z
          _h_dsigdpt_z->fill(pmom.pT()/GeV, weight);
          _h_dsigdpt_scaled_z->fill(pmom.pT()/GeV * MW_MZ, weight);
        }
      } else {
        // There is no Z -> ee candidate... so this must be a W event, right?
        const LeadingParticlesFinalState& enuFS = applyProjection<LeadingParticlesFinalState>(event, "enuFS");
        const LeadingParticlesFinalState& enubFS = applyProjection<LeadingParticlesFinalState>(event, "enubFS");
        static size_t Wcount = 0;

        // Fill W pT distributions
        ParticleVector Wdaughters;
        if (enuFS.particles().size() == 2 && enubFS.empty()) {
          Wdaughters = enuFS.particles();
        } else if (enuFS.empty() && enubFS.particles().size() == 2) {
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



    void finalize() {
      // Get cross-section per event (i.e. per unit weight) from generator
      const double xSecPerEvent = crossSectionPerEvent()/picobarn;

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
        const double BRZEE_BRWENU = 0.033632 / 0.1073;
        const double scalefactor = (xSecW / wpt_integral) / (xSecZ / zpt_scaled_integral) * BRZEE_BRWENU;
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


    //@}
 
  private:

    /// @name Event counters for cross section normalizations
    //@{
    double _eventsFilledW;
    double _eventsFilledZ;
    //@}

    //@{
    /// Histograms
    AIDA::IHistogram1D* _h_dsigdpt_w;
    AIDA::IHistogram1D* _h_dsigdpt_z;
    AIDA::IHistogram1D* _h_dsigdpt_scaled_z;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<D0_2001_S4674421> plugin_D0_2001_S4674421;

}
