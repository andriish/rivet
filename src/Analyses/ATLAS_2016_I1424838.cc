// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FParameter.hh"
#include "Rivet/Projections/Spherocity.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include <sstream>
#include <assert.h>     /* assert */
namespace Rivet {


class ATLAS_2016_I1424838 : public Analysis {
public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2016_I1424838()
        : Analysis("ATLAS_2016_I1424838")
    {


    }

    //@}


public:

    // Convenience method for histogram booking
    string mkHistoName(int idDS, int channel, int i) {
      std::stringstream s;
      s << "d0"  << idDS << "-x0" << channel  << "-y0" << i+1;
      return s.str();
    }
    
    /// Book histograms and initialise projections before the run
    void init() {

       
      
      // Charged particles inside acceptance region
      const ChargedFinalState cfs(Cuts::abseta<2.5 && Cuts::pT> 500*MeV);
      addProjection(cfs, "CFS");

      // ZFinders
      ZFinder zfinder(cfs, Cuts::abseta<2.4 && Cuts::pT>20.0*GeV, PID::ELECTRON, 66*GeV, 116*GeV, 0.1, ZFinder::CLUSTERNODECAY);
      //ZFinder zfinder(   -2.4, 2.4, 20.0, ELECTRON, 66.0*GeV, 116.0*GeV, 0.1, true, true);
      addProjection(zfinder, "ZFinder");
      ZFinder zfinder_mu(cfs, Cuts::abseta<2.4 && Cuts::pT>20.0*GeV, PID::MUON, 66*GeV, 116*GeV, 0.1, ZFinder::CLUSTERNODECAY);
      //ZFinder zfinder_mu(-2.4, 2.4, 20.0, MUON,     66.0*GeV, 116.0*GeV, 0.1, true, true);
      addProjection(zfinder_mu, "ZFinderMu");

      // This CFS only contains charged particles inside the acceptance excluding the leptons
      VetoedFinalState remfs(cfs);
      remfs.addVetoOnThisFinalState(zfinder);
      remfs.addVetoOnThisFinalState(zfinder_mu);
      addProjection(remfs, "REMFS");
      
      //const FParameter fparam(zfinder.remainingFinalState());
      const FParameter fparam(remfs);
      addProjection(fparam, "FParameter_");

      //const Spherocity sphero(zfinder.remainingFinalState());
      const Spherocity sphero(remfs);
      addProjection(sphero, "Spherocity_");


      // Booking of ES histos
      for (size_t alg = 0; alg < 5; ++alg) {
        // Book the inclusive histograms
        Elec_Ntrk[alg]         = bookHisto1D(mkHistoName(1, 1, alg));
        Elec_SumPt[alg]        = bookHisto1D(mkHistoName(2, 1, alg));
        Elec_Beamthrust[alg]   = bookHisto1D(mkHistoName(3, 1, alg)); 
        Elec_Thrust[alg]       = bookHisto1D(mkHistoName(4, 1, alg));
        Elec_FParam[alg]       = bookHisto1D(mkHistoName(5, 1, alg));  
        Elec_Spherocity[alg]   = bookHisto1D(mkHistoName(6, 1, alg));  
        Muon_Ntrk[alg]         = bookHisto1D(mkHistoName(1, 2, alg));
        Muon_SumPt[alg]        = bookHisto1D(mkHistoName(2, 2, alg));
        Muon_Beamthrust[alg]   = bookHisto1D(mkHistoName(3, 2, alg)); 
        Muon_Thrust[alg]       = bookHisto1D(mkHistoName(4, 2, alg));
        Muon_FParam[alg]       = bookHisto1D(mkHistoName(5, 2, alg));  
        Muon_Spherocity[alg]   = bookHisto1D(mkHistoName(6, 2, alg));  
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
     
        // Get generator weight
        const double weight = event.weight();

        // Check for Z boson in event
        const ZFinder& zfinder    = applyProjection<ZFinder>(event, "ZFinder");
        const ZFinder& zfinder_mu = applyProjection<ZFinder>(event, "ZFinderMu");
        bool isElec = false;
        bool isMuon = false;
        if (zfinder.bosons().size() != 1) {
            MSG_DEBUG("Num e+ e- pairs found = " << zfinder.bosons().size());
        }
        else {
          isElec = true;
        }
        
        if (zfinder_mu.bosons().size() != 1) {
            MSG_DEBUG("Num mu+ mu- pairs found = " << zfinder_mu.bosons().size());
        }
        else {
          isMuon = true;
        }

        // Only accept events with exactly two electrons or exactly two muons
        if (isElec && isMuon) vetoEvent;
        if (!(isElec || isMuon)) vetoEvent;

        // this determines the Zpt phase-space
        double zpT=-1000;
        if (isElec) zpT = zfinder.bosons()[0].momentum().perp();
        if (isMuon) zpT = zfinder_mu.bosons()[0].momentum().perp();

        unsigned int alg;
        if ( zpT/MeV < 6000.) alg=1;
        else if ( zpT/MeV >=  6000. &&  zpT/MeV < 12000.) alg=2;
        else if ( zpT/MeV >= 12000. &&  zpT/MeV < 25000.) alg=3;
        else alg=4;


        assert(alg<5);
        assert(alg>0);

        // All charged particles within |eta|<2.5 except the leptons from Z-decay
        const VetoedFinalState& remfs = applyProjection<VetoedFinalState>(event, "REMFS");
        // sumPt and Beamthrust (the latter will only be filled if the min Nch criterion is met)
        // and Thrust preparation
        
        double sumPt      = 0.0;
        double beamThrust = 0.0;

       

        std::vector<Vector3> momenta;
        foreach(const Particle& p , remfs.particles()) {
          double pT = fabs(p.momentum().pT());
          double eta = p.momentum().eta();
          sumPt += pT;
          beamThrust += pT*exp(-1.0*fabs(eta));
          Vector3 mom = p.momentum().vector3();
          mom.setZ(0.0);
          momenta.push_back(mom);
        }
     

        // Fill inclusive histos
        if (isElec) {
          Elec_Ntrk[alg]       ->fill(remfs.size(),        weight);
          Elec_Ntrk[0]         ->fill(remfs.size(),        weight);
          Elec_SumPt[alg]      ->fill(sumPt,               weight);
          Elec_SumPt[0]        ->fill(sumPt,               weight);
        }
        if (isMuon) {
          Muon_Ntrk[alg]       ->fill(remfs.size(),        weight);
          Muon_Ntrk[0]         ->fill(remfs.size(),        weight);
          Muon_SumPt[alg]      ->fill(sumPt,               weight);
          Muon_SumPt[0]        ->fill(sumPt,               weight);
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
          const FParameter& fparam = applyProjection<FParameter>(event, "FParameter_");
          // Spherocity
          const Spherocity& sphero = applyProjection<Spherocity>(event, "Spherocity_");
          
          // Histos differential in NMPI 

          // Fill inclusive histos
          if (isElec) {
            Elec_Thrust[alg]     ->fill(thrust,              weight);
            Elec_Thrust[0]       ->fill(thrust,              weight);
            Elec_FParam[alg]     ->fill(fparam.F(),          weight);
            Elec_FParam[0]       ->fill(fparam.F(),          weight);
            Elec_Spherocity[alg] ->fill(sphero.spherocity(), weight);
            Elec_Spherocity[0]   ->fill(sphero.spherocity(), weight);
            Elec_Beamthrust[alg] ->fill(beamThrust,          weight);
            Elec_Beamthrust[0]   ->fill(beamThrust,          weight);
          }
          if (isMuon) {
            Muon_Thrust[alg]     ->fill(thrust,              weight);
            Muon_Thrust[0]       ->fill(thrust,              weight);
            Muon_FParam[alg]     ->fill(fparam.F(),          weight);
            Muon_FParam[0]       ->fill(fparam.F(),          weight);
            Muon_Spherocity[alg] ->fill(sphero.spherocity(), weight);
            Muon_Spherocity[0]   ->fill(sphero.spherocity(), weight);
            Muon_Beamthrust[alg] ->fill(beamThrust,          weight);
            Muon_Beamthrust[0]   ->fill(beamThrust,          weight);
          }
        }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
        for (size_t alg = 0; alg < 5; ++alg) {
          normalize(Elec_Ntrk[alg]);
          normalize(Elec_SumPt[alg]);
          normalize(Elec_Beamthrust[alg]);
          normalize(Elec_Thrust[alg]);
          normalize(Elec_FParam[alg]);
          normalize(Elec_Spherocity[alg]);
          normalize(Muon_Ntrk[alg]);
          normalize(Muon_SumPt[alg]);
          normalize(Muon_Beamthrust[alg]);
          normalize(Muon_Thrust[alg]);
          normalize(Muon_FParam[alg]);
          normalize(Muon_Spherocity[alg]);
        }
    }

private:

     Histo1DPtr Elec_Ntrk[5];
     Histo1DPtr Elec_SumPt[5];
     Histo1DPtr Elec_Beamthrust[5];
     Histo1DPtr Elec_Thrust[5];
     Histo1DPtr Elec_FParam[5];
     Histo1DPtr Elec_Spherocity[5];

     Histo1DPtr Muon_Ntrk[5];
     Histo1DPtr Muon_SumPt[5];
     Histo1DPtr Muon_Beamthrust[5];
     Histo1DPtr Muon_Thrust[5];
     Histo1DPtr Muon_FParam[5];
     Histo1DPtr Muon_Spherocity[5];
};

// This global object acts as a hook for the plugin system
DECLARE_RIVET_PLUGIN(ATLAS_2016_I1424838);

}
