// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/EventMixingFinalState.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class EvMixing_TestAnalysis : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(EvMixing_TestAnalysis);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs(Cuts::abseta < 0.8 && Cuts::pT > 150*MeV);
      declare(fs, "FS");
      const FinalState fsmix(Cuts::abseta < 0.8 && Cuts::pT > 150*MeV);
      declare(fsmix, "MIX");

      // Declare centrality projection
      const CentralityProjection& cp = declareCentrality(ALICE::V0MMultiplicity(),
        "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Declare event mixing projection
      declare(EventMixingCentrality(&cp, fsmix, 10, 0, 100, 10), "EVMIXFS");

      // Book histograms
      for (int i = 0; i < NBINS; i++) {
        _hDEtaDPhiSame[i] = bookHisto2D("detadphisame" + toString(i), 36, -0.5*pi, 1.5*pi, 32, -1.6, 1.6,
          "Same", "$\\Delta\\phi$ (rad)", "$\\Delta\\eta$",
          "$1 / N_{same} {\\rm d}N_{same} / {\\rm d}\\Delta\\eta{\\rm d}\\Delta\\phi$ (rad$^-1$)");
        _hDEtaDPhiMixed[i] = bookHisto2D("detadphimixed" + toString(i), 36, -0.5*pi, 1.5*pi, 32, -1.6, 1.6,
          "Mixed", "$\\Delta\\phi$ (rad)", "$\\Delta\\eta$",
          "$1 / N_{mixed} {\\rm d}N_{mixed} / {\\rm d}\\Delta\\eta{\\rm d}\\Delta\\phi$ (rad$^-1$)");
        _hDEtaDPhi[i] = bookHisto2D("detadphi" + toString(i), 36, -0.5*pi, 1.5*pi, 32, -1.6, 1.6,
          "DEtaDPhi", "$\\Delta\\phi$ (rad)", "$\\Delta\\eta$",
          "$1 / N {\\rm d}N / {\\rm d}\\Delta\\eta{\\rm d}\\Delta\\phi$ (rad$^-1$)");
      }

    }
    /// @brief Calculate angular distance between particles.
    double deltaPhi(double a1, double a2){
      double dif = a1 - a2;
      while (dif < -M_PI/2)
        dif += 2*M_PI;
      while (dif > 3*M_PI/2)
        dif -= 2*M_PI;
      return dif;
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get centrality bin number for current event
      int nr = (int)(apply<CentralityProjection>(event, "V0M")() / 10);
      if (nr > NBINS) vetoEvent;

      // Get particles from current event
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      Particles particles = fs.particles();

      // Fill 'same' histograms
      for (const Particle& particle1 : particles) {
        for (const Particle& particle2 : particles) {
          if (particle1 != particle2) {
            double dPhi = deltaPhi(particle1.phi(), particle2.phi());
            if (dPhi < -0.5 * pi) dPhi += 2 * pi;
            double dEta = abs(particle1.eta()-particle2.eta());
            _hDEtaDPhiSame[nr]->fill(dPhi, dEta, event.weight());
          }
        }
      }

      // Get events for mixing
      const EventMixingCentrality& evmixfs = applyProjection<EventMixingCentrality>(event, "EVMIXFS");
      vector<Particles> mixedEvents = evmixfs.getMixingEvents();

      // Fill 'mixed' histograms
      if (mixedEvents.size() > 0) { // not really needed
        for (const auto& mixedParticles : mixedEvents) {
          //std::cout << particles.size() << "\t" << mixedParticles.size() << std::endl;
          for (const Particle& particle1 : particles) {
            for (const Particle& particle2 : mixedParticles) {
              if (particle1 != particle2) {
                double dPhi = deltaPhi(particle1.phi(), particle2.phi());
                if (dPhi < -0.5 * pi) dPhi += 2 * pi;
                double dEta = abs(particle1.eta()-particle2.eta());
                _hDEtaDPhiMixed[nr]->fill(dPhi, dEta, 1. / mixedEvents.size());
              }
            }
          }
        }
        //std::cout << std::endl;
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for (int ibin = 0; ibin < NBINS; ibin++) {
        if (_hDEtaDPhiMixed[ibin]->numEntries() < 1)
          continue;
        std::cout << "Creating final histogram " << ibin << std::endl;
        // We have to do that manually (bin by bin) because there is no way to use divide() method on Histo2DPtr
        // and it is not possible to use Scatter3DPtr because there is no bookScatter3DPtr method
        for (unsigned int ixy = 0; ixy < _hDEtaDPhiSame[ibin]->numBins(); ixy++) {
          _hDEtaDPhi[ibin]->fillBin(ixy, _hDEtaDPhiSame[ibin]->bin(ixy).sumW() / _hDEtaDPhiMixed[ibin]->bin(ixy).sumW());
        }
      }

    }

    //@}

  private:

    static const int NBINS = 10;

    /// @name Histograms
    //@{
    Histo2DPtr _hDEtaDPhiSame[NBINS];
    Histo2DPtr _hDEtaDPhiMixed[NBINS];
    Histo2DPtr _hDEtaDPhi[NBINS];
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EvMixing_TestAnalysis);


}
