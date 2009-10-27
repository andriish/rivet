// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  /// @brief Measurement differential Z/gamma* + jet +X cross sections
  /// @author Gavin Hesketh, Andy Buckley, Frank Siegert
  class D0_2008_S7863608 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    D0_2008_S7863608() : Analysis("D0_2008_S7863608")
    {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
    }
    
    //@}


    /// @name Analysis methods
    //@{     
    
    /// Book histograms
    void init() {
      ZFinder zfinder(-1.7, 1.7, 15.0*GeV, ELECTRON, 65.0*GeV, 115.0*GeV, 0.2);
      addProjection(zfinder, "ZFinder");
      
      FastJets conefinder(zfinder.remainingFinalState(), FastJets::D0ILCONE, 0.5, 20.0*GeV);
      addProjection(conefinder, "ConeFinder");

      _h_jet_pT_cross_section = bookHistogram1D(1, 1, 1);
      _h_jet_y_cross_section = bookHistogram1D(2, 1, 1);
      _h_Z_pT_cross_section = bookHistogram1D(3, 1, 1);
      _h_Z_y_cross_section = bookHistogram1D(4, 1, 1);
      _h_total_cross_section = bookHistogram1D(5, 1, 1);  
    }
    
    

    // Do the analysis 
    void analyze(const Event& e) {
      const double weight = e.weight();
      
      const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
      if (zfinder.particles().size()==1) {
        const JetAlg& jetpro = applyProjection<JetAlg>(e, "ConeFinder");
        const Jets& jets = jetpro.jetsByPt(20.0*GeV);
        Jets jets_cut;
        foreach (const Jet& j, jets) {
          if (fabs(j.momentum().pseudorapidity()) < 2.8) {
            jets_cut.push_back(j);
          }
        }
        
        // Return if there are no jets:
        if(jets_cut.size()<1) {
          getLog() << Log::DEBUG << "Skipping event " << e.genEvent().event_number()
                   << " because no jets pass cuts " << endl;
          vetoEvent;
        }
        
        // cut on Delta R between jet and muons
        foreach (const Jet& j, jets_cut) {
          foreach (const Particle& mu, zfinder.constituentsFinalState().particles()) {
            if (deltaR(mu.momentum().pseudorapidity(), mu.momentum().azimuthalAngle(),
                       j.momentum().pseudorapidity(), j.momentum().azimuthalAngle()) < 0.5) {
              vetoEvent;
            }
          }
        }
        
        const FourMomentum Zmom = zfinder.particles()[0].momentum();
        
        // In jet pT
        _h_jet_pT_cross_section->fill( jets_cut[0].momentum().pT(), weight);
        _h_jet_y_cross_section->fill( fabs(jets_cut[0].momentum().rapidity()), weight);
        
        // In Z pT
        _h_Z_pT_cross_section->fill(Zmom.pT(), weight);
        _h_Z_y_cross_section->fill(fabs(Zmom.rapidity()), weight);
        
        _h_total_cross_section->fill(1960.0, weight);
      }
    }
    
    
    
    /// Finalize
    void finalize() {
      const double invlumi = crossSection()/sumOfWeights();
      scale(_h_total_cross_section, invlumi);
      scale(_h_jet_pT_cross_section, invlumi);
      scale(_h_jet_y_cross_section, invlumi);
      scale(_h_Z_pT_cross_section, invlumi);
      scale(_h_Z_y_cross_section, invlumi);
    }
    
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_jet_pT_cross_section;
    AIDA::IHistogram1D * _h_jet_y_cross_section;
    AIDA::IHistogram1D * _h_Z_pT_cross_section;
    AIDA::IHistogram1D * _h_Z_y_cross_section;
    AIDA::IHistogram1D * _h_total_cross_section;
    //@}

  };

    
    
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<D0_2008_S7863608> plugin_D0_2008_S7863608;
  
}
