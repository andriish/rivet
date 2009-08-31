// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  /// @brief Measurement of isolated gamma + jet + X differential cross-sections
  /// Inclusive isolated gamma + jet cross-sections, differential in pT(gamma), for 
  /// various photon and jet rapidity bins.
  ///
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7719523 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_2008_S7719523()
      : Analysis("D0_2008_S7719523")
    {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
      
      // General FS
      FinalState fs(-5.0, 5.0);
      addProjection(fs, "FS");
      
      // Get leading photon
      LeadingParticlesFinalState photonfs(fs, -1.0, 1.0);
      photonfs.addParticleId(PHOTON);
      addProjection(photonfs, "LeadingPhoton");
      
      // FS for jets excludes the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      addProjection(vfs, "JetFS");
    } 
    
    //@}


    /// @name Analysis methods
    //@{ 
    
    /// Book histograms
    void init() {
      _h_central_same_cross_section = bookHistogram1D(1, 1, 1);
      _h_central_opp_cross_section  = bookHistogram1D(2, 1, 1);
      _h_forward_same_cross_section = bookHistogram1D(3, 1, 1);
      _h_forward_opp_cross_section  = bookHistogram1D(4, 1, 1); 
    }
    
    

    /// Do the analysis 
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get the photon
      const FinalState& photonfs = applyProjection<FinalState>(event, "LeadingPhoton");
      if (photonfs.particles().size() != 1) {
        getLog() << Log::DEBUG << "No photon found" << endl;
        vetoEvent;
      }
      const FourMomentum photon = photonfs.particles().front().momentum();
      if (photon.pT()/GeV < 30) {
        getLog() << Log::DEBUG << "Leading photon has pT < 30 GeV: " << photon.pT()/GeV << endl;
        vetoEvent;
      }
      
      // Get all charged particles
      const FinalState& fs = applyProjection<FinalState>(event, "JetFS");
      if (fs.isEmpty()) {
        vetoEvent;
      }
      
      // Isolate photon by ensuring that a 0.4 cone around it contains less than 7% of the photon's energy
      const double egamma = photon.E();
      double econe = 0.0;
      foreach (const Particle& p, fs.particles()) {
        if (deltaR(photon, p.momentum()) < 0.4) {
          econe += p.momentum().E();
          // Veto as soon as E_cone gets larger
          if (econe/egamma > 0.07) {
            getLog() << Log::DEBUG << "Vetoing event because photon is insufficiently isolated" << endl;
            vetoEvent;
          }
        }
      }
      
      
      /// @todo Allow proj creation w/o FS as ctor arg, so that calc can be used more easily.
      FastJets jetpro(fs, FastJets::D0ILCONE, 0.7); //< @todo This fs arg makes no sense!
      jetpro.calc(fs.particles());
      Jets isolated_jets;
      foreach (const Jet& j, jetpro.jets()) {
        const FourMomentum pjet = j.momentum();
        const double dr = deltaR(photon.pseudorapidity(), photon.azimuthalAngle(),
                                 pjet.pseudorapidity(), pjet.azimuthalAngle());
        if (dr > 0.7 && pjet.pT()/GeV > 15) {
          isolated_jets.push_back(j);
        }
      }
      
      getLog() << Log::DEBUG << "Num jets after isolation and pT cuts = " 
               << isolated_jets.size() << endl;
      if (isolated_jets.empty()) {
        getLog() << Log::DEBUG << "No jets pass cuts" << endl;
        vetoEvent;
      }
      
      // Sort by pT and get leading jet
      sort(isolated_jets.begin(), isolated_jets.end(), cmpJetsByPt);
      const FourMomentum leadingJet = isolated_jets.front().momentum();
      int photon_jet_sign = sign( leadingJet.rapidity() * photon.rapidity() );
      
      // Veto if leading jet is outside plotted rapidity regions
      const double abs_y1 = fabs(leadingJet.rapidity());
      if (inRange(abs_y1, 0.8, 1.5) || abs_y1 > 2.5) {
        getLog() << Log::DEBUG << "Leading jet falls outside acceptance range; |y1| = " 
                 << abs_y1 << endl;
        vetoEvent;
      }
      
      // Fill histos
      if (fabs(leadingJet.rapidity()) < 0.8) { 
        if (photon_jet_sign >= 1) {
          _h_central_same_cross_section->fill(photon.pT(), weight);
        } else {
          _h_central_opp_cross_section->fill(photon.pT(), weight);
        }
      } else if (inRange( fabs(leadingJet.rapidity()), 1.5, 2.5)) {
        if (photon_jet_sign >= 1) {
          _h_forward_same_cross_section->fill(photon.pT(), weight);
        } else {
          _h_forward_opp_cross_section->fill(photon.pT(), weight); 
        }
      }
      
    }
    
    
    
    /// Finalize
    void finalize() {
      const double lumi_gen = sumOfWeights()/crossSection();
      const double dy_photon = 2.0;
      const double dy_jet_central = 1.6;
      const double dy_jet_forward = 2.0;
      
      // Cross-section ratios (6 plots)
      // Central/central and forward/forward ratios
      AIDA::IHistogramFactory& hf = histogramFactory();
      const string dir = histoDir();
      
      hf.divide(dir + "/d05-x01-y01", *_h_central_opp_cross_section, *_h_central_same_cross_section);
      hf.divide(dir + "/d08-x01-y01", *_h_forward_opp_cross_section, *_h_forward_same_cross_section);
      
      // Central/forward ratio combinations
      hf.divide(dir + "/d06-x01-y01", *_h_central_same_cross_section,
                *_h_forward_same_cross_section)->scale(dy_jet_forward/dy_jet_central, 1);
      hf.divide(dir + "/d07-x01-y01", *_h_central_opp_cross_section,
                *_h_forward_same_cross_section)->scale(dy_jet_forward/dy_jet_central, 1);
      hf.divide(dir + "/d09-x01-y01", *_h_central_same_cross_section,
                *_h_forward_opp_cross_section)->scale(dy_jet_forward/dy_jet_central, 1);
      hf.divide(dir + "/d10-x01-y01", *_h_central_opp_cross_section,
                *_h_forward_opp_cross_section)->scale(dy_jet_forward/dy_jet_central, 1);
      
      // Use generator cross section for remaining histograms
      scale(_h_central_same_cross_section, 1.0/lumi_gen * 1.0/dy_photon * 1.0/dy_jet_central);
      scale(_h_central_opp_cross_section, 1.0/lumi_gen * 1.0/dy_photon * 1.0/dy_jet_central);
      scale(_h_forward_same_cross_section, 1.0/lumi_gen * 1.0/dy_photon * 1.0/dy_jet_forward);
      scale(_h_forward_opp_cross_section, 1.0/lumi_gen * 1.0/dy_photon * 1.0/dy_jet_forward);
    }
    
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _h_central_same_cross_section;
    AIDA::IHistogram1D* _h_central_opp_cross_section;
    AIDA::IHistogram1D* _h_forward_same_cross_section;
    AIDA::IHistogram1D* _h_forward_opp_cross_section;
    //@}

  };

    
    
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<D0_2008_S7719523> plugin_D0_2008_S7719523;
  
}
