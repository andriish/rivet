// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

#include "Rivet/Particle.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "HepMC/GenEvent.h"


namespace Rivet {


  /// @brief ATLAS inclusive b-jet pT spectrum, di-jet Mass and di-jet chi
  class ATLAS_2011_I930220: public Analysis {
  public:

    ATLAS_2011_I930220() : Analysis("ATLAS_2011_I930220")
    {    }


    void init() {
      FinalState fs(-3.5, 3.5);
      addProjection(fs, "FinalState");

      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      fj.useInvisibles();
      addProjection(fj, "Jets");

      double ybins[] = {0.0, 0.3, 0.8, 1.2, 2.1};
      for (size_t i = 0; i < 4; ++i) {
        _bjetpT_SV0.addHistogram(ybins[i], ybins[i+1], bookHistogram1D(i+1, 1, 1));
      }

      _bjetpT_SV0_All    = bookHistogram1D(5, 1, 1);
      _bjetpT_pTRel      = bookHistogram1D(6, 1, 1);
      _dijet_mass        = bookHistogram1D(7, 1, 1);
      _dijet_phi         = bookHistogram1D(8, 1, 1);
      _dijet_chi_110_370 = bookHistogram1D(9, 1, 1);
      _dijet_chi_370_850 = bookHistogram1D(10, 1, 1);

      _chiCounter1 = 0.0;
      _chiCounter2 = 0.0;
      _phiCounter  = 0.0;
    }


    void analyze(const Event& evt) {
      const double weight = evt.weight();

      // Finding all b-hadrons within event
      /// @todo Rewrite via UnstableFinalState or HeavyHadrons
      std::vector<HepMC::GenParticle*> bHadrons;
      std::vector<HepMC::GenParticle*> allParticles = particles(evt.genEvent());
      foreach (HepMC::GenParticle* p, allParticles) {
        if ( !Rivet::PID::isHadron(p->pdg_id()) || !Rivet::PID::hasBottom(p->pdg_id()) ) continue;
        if ( p->momentum().perp() < 5*GeV ) continue;
        bHadrons.push_back(p);
      }

      const Jets jets = applyProjection<JetAlg>(evt, "Jets").jetsByPt(15.0*GeV);

      FourMomentum leadingJet, subleadingJet;
      int leadJet = 0, subJet = 0;
      foreach (const Jet& j, jets) {
        const FourMomentum jet = j.momentum();
        bool hasB = false;
        /// @todo Convert this to Rivet Particle
        foreach (HepMC::GenParticle* bgp, bHadrons) {
          HepMC::FourVector b = bgp->momentum();
          double dR = deltaR(jet.eta(), jet.phi(), b.eta(), b.phi());
          if (dR < 0.3) {
            hasB = true;
            break;
          }
        }

        // Identify and classify the leading and subleading jets
        if (fabs(jet.rapidity()) < 2.1) { //< Move this into the jets defn
          if (!leadJet) {
            leadingJet = jet;
            leadJet = (hasB && jet.pT() > 40*GeV) ? 2 : 1;
          }
          if (leadJet && !subJet) {
            subleadingJet = jet;
            subJet = (hasB && jet.pT() > 40*GeV) ? 2 : 1;
          }
          if (hasB) {
            _bjetpT_SV0.fill(fabs(jet.rapidity()), jet.pT()/GeV, weight);
            _bjetpT_SV0_All->fill(jet.pT()/GeV, weight);
            _bjetpT_pTRel->fill(jet.pT()/GeV, weight);
          }
        }
      }

      // Di-b-jet plots require both the leading and subleading jets to be b-tagged and have pT > 40 GeV
      if (leadJet == 2 && subJet == 2) {
        const double mass = FourMomentum( leadingJet + subleadingJet ).mass();
        _dijet_mass->fill(mass/GeV, weight);

        // Plot dphi for high-mass di-b-jets
        if (mass > 110*GeV) {
          _phiCounter += weight;
          const double d_phi = deltaPhi( leadingJet.phi(), subleadingJet.phi() );
          _dijet_phi->fill(fabs(d_phi), weight);
        }

        // Plot chi for low y_boost di-b-jets (in two high-mass bins)
        const double y_boost = 0.5 * (leadingJet.rapidity() + subleadingJet.rapidity());
        const double chi = exp( fabs( leadingJet.rapidity() - subleadingJet.rapidity() ) );
        if ( fabs(y_boost) < 1.1 ) {
          if (inRange(mass/GeV, 110, 370)) {
            _chiCounter1 += weight;
            _dijet_chi_110_370->fill(chi, weight);
          } else if (inRange(mass/GeV, 370, 850)) {
            _chiCounter2 += weight;
            _dijet_chi_370_850->fill(chi, weight);
          }
        }
      }
    }


    void finalize() {
      // Normalizing to cross-section and mass
      // Additional factors represent the division by rapidity
      const double xsec = crossSectionPerEvent()/(picobarn);
      const double chiScale1 = 1 / _chiCounter1 / 260.0;
      const double chiScale2 = 1 / _chiCounter2 / 480.0;
      const double phiScale  = 1 / _phiCounter;

      _bjetpT_SV0.scale(xsec/2, this);
      scale(_bjetpT_SV0_All, xsec);
      scale(_bjetpT_pTRel, xsec);
      scale(_dijet_mass, xsec);
      scale(_dijet_phi, phiScale );
      scale(_dijet_chi_110_370, chiScale1);
      scale(_dijet_chi_370_850, chiScale2);
    }


  private:

    BinnedHistogram<double> _bjetpT_SV0;

    AIDA::IHistogram1D* _bjetpT_SV0_All;
    AIDA::IHistogram1D* _bjetpT_pTRel;
    AIDA::IHistogram1D* _dijet_mass;
    AIDA::IHistogram1D* _dijet_phi;
    AIDA::IHistogram1D* _dijet_chi_110_370;
    AIDA::IHistogram1D* _dijet_chi_370_850;

    double _chiCounter1;
    double _chiCounter2;
    double _phiCounter;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I930220);

}
