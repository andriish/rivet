#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FParameter.hh"
#include "Rivet/Projections/Spherocity.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Event shapes in leptonic $Z$-events
  class ATLAS_2016_I1424838 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1424838);


    /// Book histograms and initialise projections before the run
    void init() {

      // Charged particles inside acceptance region
      const ChargedFinalState cfs(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV);
      declare(cfs, "CFS");

      // ZFinders
      ZFinder zfinder(cfs, Cuts::abseta<2.4 && Cuts::pT>20.0*GeV, PID::ELECTRON, 66*GeV, 116*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder, "ZFinder");
      ZFinder zfinder_mu(cfs, Cuts::abseta<2.4 && Cuts::pT>20.0*GeV, PID::MUON, 66*GeV, 116*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder_mu, "ZFinderMu");

      // This CFS only contains charged particles inside the acceptance excluding the leptons
      VetoedFinalState remfs(cfs);
      remfs.addVetoOnThisFinalState(zfinder);
      remfs.addVetoOnThisFinalState(zfinder_mu);
      declare(remfs, "REMFS");

      const FParameter fparam(remfs);
      declare(fparam, "FParameter_");

      const Spherocity sphero(remfs);
      declare(sphero, "Spherocity_");


      // Booking of ES histos
      for (size_t alg = 0; alg < 5; ++alg) {
        // Book the inclusive histograms
        size_t offset = alg * 6;
        book(_h_Elec_Ntrk[alg]         , offset+1,  1, 1);
        book(_h_Elec_SumPt[alg]        , offset+2,  1, 1);
        book(_h_Elec_Beamthrust[alg]   , offset+3,  1, 1);
        book(_h_Elec_Thrust[alg]       , offset+4,  1, 1);
        book(_h_Elec_FParam[alg]       , offset+5,  1, 1);
        book(_h_Elec_Spherocity[alg]   , offset+6,  1, 1);
        book(_h_Muon_Ntrk[alg]         , offset+31, 1, 1);
        book(_h_Muon_SumPt[alg]        , offset+32, 1, 1);
        book(_h_Muon_Beamthrust[alg]   , offset+33, 1, 1);
        book(_h_Muon_Thrust[alg]       , offset+34, 1, 1);
        book(_h_Muon_FParam[alg]       , offset+35, 1, 1);
        book(_h_Muon_Spherocity[alg]   , offset+36, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Check for Z boson in event
      const ZFinder& zfinder    = apply<ZFinder>(event, "ZFinder");
      MSG_DEBUG("Num e+ e- pairs found = " << zfinder.bosons().size());
      const bool isElec = zfinder.bosons().size() == 1;

      const ZFinder& zfinder_mu = apply<ZFinder>(event, "ZFinderMu");
      MSG_DEBUG("Num mu+ mu- pairs found = " << zfinder_mu.bosons().size());
      const bool isMuon = zfinder_mu.bosons().size() == 1;

      // Only accept events with exactly two electrons or exactly two muons
      if (isElec && isMuon) vetoEvent;
      if (!(isElec || isMuon)) vetoEvent;

      // This determines the Zpt phase-space
      double zpT = -1000;
      if (isElec) zpT = zfinder.bosons()[0].pT();
      if (isMuon) zpT = zfinder_mu.bosons()[0].pT();

      unsigned int alg = 4; //< for > 25 GeV
      if (zpT < 6*GeV) alg = 1;
      else if (inRange(zpT/GeV, 6, 12)) alg = 2;
      else if (inRange(zpT/GeV, 12, 25)) alg = 3;
      assert(alg < 5);
      assert(alg > 0);

      // All charged particles within |eta|<2.5 except the leptons from Z-decay
      const VetoedFinalState& remfs = apply<VetoedFinalState>(event, "REMFS");
      // sumPt and Beamthrust (the latter will only be filled if the min Nch criterion is met)
      // and Thrust preparation
      double sumPt = 0.0, beamThrust = 0.0;
      vector<Vector3> momenta;
      for (const Particle& p : remfs.particles()) {
        const double pT = p.pT();
        sumPt += pT;
        beamThrust += pT*exp(-p.abseta());
        const Vector3 mom = p.mom().pTvec();
        momenta.push_back(mom);
      }

      // Fill inclusive histos
      if (isElec) {
        _h_Elec_Ntrk[alg]       ->fill(remfs.size());
        _h_Elec_Ntrk[0]         ->fill(remfs.size());
        _h_Elec_SumPt[alg]      ->fill(sumPt);
        _h_Elec_SumPt[0]        ->fill(sumPt);
      }
      if (isMuon) {
        _h_Muon_Ntrk[alg]       ->fill(remfs.size());
        _h_Muon_Ntrk[0]         ->fill(remfs.size());
        _h_Muon_SumPt[alg]      ->fill(sumPt);
        _h_Muon_SumPt[0]        ->fill(sumPt);
      }

      // Skip event shape calculation if we don't match the minimum Nch criterion
      if (remfs.size() >=2) {

        // Eventshape calculations

        // Calculate transverse Thrust using all charged FS particles except the lepton
        // This is copied/inspired from the CMS_6000011_S8957746 analysis
        if (momenta.size() == 2) {
          // We need to use a ghost so that Thrust.calc() doesn't return 1.
          momenta.push_back(Vector3(1e-10*MeV, 0., 0.));
        }
        Thrust thrustC;
        thrustC.calc(momenta);

        double thrust = thrustC.thrust();

        // F-Parameter
        const FParameter& fparam = apply<FParameter>(event, "FParameter_");
        // Spherocity
        const Spherocity& sphero = apply<Spherocity>(event, "Spherocity_");

        // Histos differential in NMPI

        // Fill inclusive histos
        if (isElec) {
          _h_Elec_Thrust[alg]     ->fill(thrust);
          _h_Elec_Thrust[0]       ->fill(thrust);
          _h_Elec_FParam[alg]     ->fill(fparam.F());
          _h_Elec_FParam[0]       ->fill(fparam.F());
          _h_Elec_Spherocity[alg] ->fill(sphero.spherocity());
          _h_Elec_Spherocity[0]   ->fill(sphero.spherocity());
          _h_Elec_Beamthrust[alg] ->fill(beamThrust/GeV);
          _h_Elec_Beamthrust[0]   ->fill(beamThrust/GeV);
        }
        if (isMuon) {
          _h_Muon_Thrust[alg]     ->fill(thrust);
          _h_Muon_Thrust[0]       ->fill(thrust);
          _h_Muon_FParam[alg]     ->fill(fparam.F());
          _h_Muon_FParam[0]       ->fill(fparam.F());
          _h_Muon_Spherocity[alg] ->fill(sphero.spherocity());
          _h_Muon_Spherocity[0]   ->fill(sphero.spherocity());
          _h_Muon_Beamthrust[alg] ->fill(beamThrust/GeV);
          _h_Muon_Beamthrust[0]   ->fill(beamThrust/GeV);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t alg = 0; alg < 5; ++alg) {
        normalize(_h_Elec_Ntrk[alg],100.);
        normalize(_h_Elec_SumPt[alg],100.);
        normalize(_h_Elec_Beamthrust[alg],100.);
        normalize(_h_Elec_Thrust[alg]);
        normalize(_h_Elec_FParam[alg]);
        normalize(_h_Elec_Spherocity[alg]);
        normalize(_h_Muon_Ntrk[alg],100.);
        normalize(_h_Muon_SumPt[alg],100.);
        normalize(_h_Muon_Beamthrust[alg],100.);
        normalize(_h_Muon_Thrust[alg]);
        normalize(_h_Muon_FParam[alg]);
        normalize(_h_Muon_Spherocity[alg]);
      }
    }


  private:

    Histo1DPtr _h_Elec_Ntrk[5];
    Histo1DPtr _h_Elec_SumPt[5];
    Histo1DPtr _h_Elec_Beamthrust[5];
    Histo1DPtr _h_Elec_Thrust[5];
    Histo1DPtr _h_Elec_FParam[5];
    Histo1DPtr _h_Elec_Spherocity[5];

    Histo1DPtr _h_Muon_Ntrk[5];
    Histo1DPtr _h_Muon_SumPt[5];
    Histo1DPtr _h_Muon_Beamthrust[5];
    Histo1DPtr _h_Muon_Thrust[5];
    Histo1DPtr _h_Muon_FParam[5];
    Histo1DPtr _h_Muon_Spherocity[5];

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1424838);


}
