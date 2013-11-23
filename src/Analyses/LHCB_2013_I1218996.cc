// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"

namespace Rivet {


  class LHCB_2013_I1218996 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    LHCB_2013_I1218996()
      : Analysis("LHCB_2013_I1218996")
    {    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections
      addProjection(UnstableFinalState(), "UFS");

      /// Book histograms
      _h_pdg411_Dplus_pT_y.addHistogram(  2.0, 2.5, bookHistogram1D(3, 1, 1));
      _h_pdg411_Dplus_pT_y.addHistogram(  2.5, 3.0, bookHistogram1D(3, 1, 2));
      _h_pdg411_Dplus_pT_y.addHistogram(  3.0, 3.5, bookHistogram1D(3, 1, 3));
      _h_pdg411_Dplus_pT_y.addHistogram(  3.5, 4.0, bookHistogram1D(3, 1, 4));
      _h_pdg411_Dplus_pT_y.addHistogram(  4.0, 4.5, bookHistogram1D(3, 1, 5));

      _h_pdg421_Dzero_pT_y.addHistogram(  2.0, 2.5, bookHistogram1D(2, 1, 1));
      _h_pdg421_Dzero_pT_y.addHistogram(  2.5, 3.0, bookHistogram1D(2, 1, 2));
      _h_pdg421_Dzero_pT_y.addHistogram(  3.0, 3.5, bookHistogram1D(2, 1, 3));
      _h_pdg421_Dzero_pT_y.addHistogram(  3.5, 4.0, bookHistogram1D(2, 1, 4));
      _h_pdg421_Dzero_pT_y.addHistogram(  4.0, 4.5, bookHistogram1D(2, 1, 5));

      _h_pdg431_Dsplus_pT_y.addHistogram(  2.0, 2.5, bookHistogram1D(5, 1, 1));
      _h_pdg431_Dsplus_pT_y.addHistogram(  2.5, 3.0, bookHistogram1D(5, 1, 2));
      _h_pdg431_Dsplus_pT_y.addHistogram(  3.0, 3.5, bookHistogram1D(5, 1, 3));
      _h_pdg431_Dsplus_pT_y.addHistogram(  3.5, 4.0, bookHistogram1D(5, 1, 4));
      _h_pdg431_Dsplus_pT_y.addHistogram(  4.0, 4.5, bookHistogram1D(5, 1, 5));

      _h_pdg413_Dstarplus_pT_y.addHistogram(  2.0, 2.5, bookHistogram1D(4, 1, 1));
      _h_pdg413_Dstarplus_pT_y.addHistogram(  2.5, 3.0, bookHistogram1D(4, 1, 2));
      _h_pdg413_Dstarplus_pT_y.addHistogram(  3.0, 3.5, bookHistogram1D(4, 1, 3));
      _h_pdg413_Dstarplus_pT_y.addHistogram(  3.5, 4.0, bookHistogram1D(4, 1, 4));
      _h_pdg413_Dstarplus_pT_y.addHistogram(  4.0, 4.5, bookHistogram1D(4, 1, 5));

      _h_pdg4122_Lambdac_pT = bookHistogram1D(1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      /// @todo Use PrimaryHadrons to avoid double counting?
      const UnstableFinalState &ufs = applyProjection<UnstableFinalState> (event, "UFS");

      foreach (const Particle& p, ufs.particles() ) {
        const PdgId pid = abs(p.pdgId());

        // We're only interested in charm hadrons
        if (!(PID::isHadron(pid) && PID::hasCharm(pid))) continue;

        // Kinematic acceptance
        /// @todo Use direct Particle methods
        const double y = fabs(p.momentum().rapidity()); //< Double analysis efficiency with a "two-sided LHCb"
        const double pT = p.momentum().pT();

        // Fiducial acceptance of the measurements
        if (pT > 8.0*GeV || y < 2.0 || y > 4.5) continue;

        /// Experimental selection removes non-prompt charm hadrons: we ignore those from b decays
        if (hasBottomAncestor(p)) continue;

        switch (pid) {
        case 411:
          _h_pdg411_Dplus_pT_y.fill(y, pT/GeV, weight);
          break;
        case 421:
          _h_pdg421_Dzero_pT_y.fill(y, pT/GeV, weight);
          break;
        case 431:
          _h_pdg431_Dsplus_pT_y.fill(y, pT/GeV, weight);
          break;
        case 413:
          _h_pdg413_Dstarplus_pT_y.fill(y, pT/GeV, weight);
          break;
        case 4122:
          _h_pdg4122_Lambdac_pT->fill(pT/GeV, weight);
          break;
        }
      }
    }


    /// Search ancestry for a b-hadron.
    ///
    /// Adapted from Rivet::Particle::hasAncestor
    ///
    /// @todo Replace with p.fromBottom() in Rivet 2.0.1
    bool hasBottomAncestor(const Particle& p) const {
      HepMC::GenVertex* prodVtx = p.genParticle().production_vertex();
      if (!prodVtx) return false;
      foreach (const HepMC::GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
        if (PID::hasBottom(ancestor->pdg_id())) return true;
      }
      return false;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double scale_factor = crossSection()/microbarn / sumOfWeights();

      /// Avoid the implicit division by the bin width in the BinnedHistogram::scale method.
      foreach(AIDA::IHistogram1D* hist, _h_pdg411_Dplus_pT_y.getHistograms()) hist->scale(scale_factor);
      foreach(AIDA::IHistogram1D* hist, _h_pdg421_Dzero_pT_y.getHistograms()) hist->scale(scale_factor);
      foreach(AIDA::IHistogram1D* hist, _h_pdg431_Dsplus_pT_y.getHistograms()) hist->scale(scale_factor);
      foreach(AIDA::IHistogram1D* hist, _h_pdg413_Dstarplus_pT_y.getHistograms()) hist->scale(scale_factor);
      _h_pdg4122_Lambdac_pT->scale(scale_factor);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    BinnedHistogram<double> _h_pdg411_Dplus_pT_y;
    BinnedHistogram<double> _h_pdg421_Dzero_pT_y;
    BinnedHistogram<double> _h_pdg431_Dsplus_pT_y;
    BinnedHistogram<double> _h_pdg413_Dstarplus_pT_y;
    AIDA::IHistogram1D *_h_pdg4122_Lambdac_pT;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCB_2013_I1218996);

}
