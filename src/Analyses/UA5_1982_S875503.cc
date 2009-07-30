#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {
  
  class UA5_1982_S875503 : public Analysis {
    
  public:
    
    /// Default constructor
    UA5_1982_S875503() : Analysis("UA5_1982_S875503") {
      const FinalState fs;
      const ChargedFinalState cfs;
      addProjection(Beam(), "Beam");
      addProjection(fs, "FS");
      addProjection(cfs, "CFS");
    }
    
    
    /// Factory method
    static Analysis* create() { 
      return new UA5_1982_S875503(); 
    }
    
    
    /// Return the name of this analysis
    string name() const {
      return "UA5_1982_S875503";
    }
    
    /// Get the SPIRES ID code
    string spiresId() const {
      return "S875503";
    }
    
    /// Get a description of the analysis.
    string description() const {
      return "";
    }
    
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "UA5";
    }
    
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1982";
    }
    
    /// Summary of analysis
    string summary() const{
      return "";
    }
    
    /// Beam conditions for this analysis
    string runInfo() const{
      return "Either PP or PPBar";
    }
    
    string collider() const{
      return "Any";
    }
    
    vector<string> authors() const{
      return vector<string>();
    }
    
    vector<string> references() const{
      return vector<string>();
    }
    
    /// @name Analysis methods
    //@{
    void init() 
    { 
    _hist_etapp = bookHistogram1D("d01-x01-y01", 35, 0., 3.5);

    _hist_etappbar = bookHistogram1D("d01-x01-y02", 35, 0., 3.5);

    _hist_nchpp = bookHistogram1D("d02-x01-y01", 18, 0., 36.); 

    _hist_nchppbar = bookHistogram1D("d02-x01-y02", 18, 0., 36.); 
    }
    
    void analyze(const Event& event) 
    {
      const Beam b = applyProjection<Beam>(event, "Beam");
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      const double weight = event.weight();
      foreach (const Particle& p, fs.particles())
      {
        if ( b.beams().first.pdgId() == b.beams().second.pdgId())
        {
          _hist_etapp->fill(p.momentum().pseudorapidity(), weight);
        }
        else if ( b.beams().first.pdgId() != b.beams().second.pdgId())
        {
          _hist_etappbar->fill(p.momentum().pseudorapidity(), weight);
        }
      }
        const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
        foreach (const Particle& p, cfs.particles())
        {
        if ( b.beams().first.pdgId() == b.beams().second.pdgId())
        {
          _hist_nchpp->fill(cfs.particles().size(), weight);
        }
        else if ( b.beams().first.pdgId() != b.beams().second.pdgId())
        {
          _hist_nchppbar->fill(cfs.particles().size(), weight);
        }
      }
    }
    
    void finalize() {
      normalize(_hist_etapp);
      normalize(_hist_etappbar);
      normalize(_hist_nchpp);
      normalize(_hist_nchppbar);      
    }
    //@}
    
  
    private:
    
    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _hist_etapp;
    AIDA::IHistogram1D* _hist_etappbar;
    AIDA::IHistogram1D* _hist_nchpp;
    AIDA::IHistogram1D* _hist_nchppbar;
    //@}


    
  };
  
//   extern "C" {
//     AnalysisBuilders getAnalysisBuilders() {
//       AnalysisBuilders fns;
//       fns["UA5_1982_S875503"] = Rivet::UA5_1982_S875503::create;
//       return fns;
//     }
//   }

  
}
