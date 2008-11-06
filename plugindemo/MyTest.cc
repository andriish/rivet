// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  class MyTest : public Analysis {
  public:

    /// Default constructor
    MyTest() {
      FinalState fs;
      ChargedFinalState cfs;
      addProjection(fs, "FS");
      addProjection(ChargedFinalState(fs), "CFS");
      //addProjection(Multiplicity(fs), "Mult");
      //addProjection(Multiplicity(cfs), "CMult");
      //addProjection(Thrust(fs), "Thrust");
    }
    
    
    /// Factory method
    static Analysis* create() { 
      return new MyTest(); 
    }
    
    
    /// Get the name of this analysis.
    string getName() const {
      return "MyTest";
    }


    // Book histograms
    void init() {
      _histTot   = bookHistogram1D("TotalMult", "Total multiplicity", 250, -0.5, 499.5);
      _histChTot = bookHistogram1D("TotalChMult", "Total charged multiplicity", 250, -0.5, 499.5);
    }
    

    // Do the analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Fill charged+neutral particle histograms
      const ParticleVector fsparts = 
        applyProjection<FinalState>(event, "FS").particles();
      unsigned int lepmult = 0;
      foreach (const Particle& p, fsparts) {
        if (PID::isLepton(p.getPdgId())) lepmult += 1;
      }
      _histTot->fill(fsparts.size(), weight);
      _histLeps->fill(lepmult, weight);
      _histNotLeps->fill(fsparts.size() - lepmult, weight);

      // Fill charged particle histograms
      const ParticleVector cfsparts = 
        applyProjection<FinalState>(event, "CFS").particles();
      unsigned int chlepmult = 0;
      foreach (const Particle& p, cfsparts) {
        if (PID::isLepton(p.getPdgId())) chlepmult += 1;
      }
      _histChTot->fill(cfsparts.size(), weight);
      _histChLeps->fill(chlepmult, weight);
      _histChNotLeps->fill(cfsparts.size() - chlepmult, weight);      
    }
    
    
    // Finalize
    void finalize() { 
      // normalize(_histTot);
      // normalize(_histChTot);
    }
    
    
  private:
    
    //@{
    /// Histograms
    AIDA::IHistogram1D *_histTot, *_histChTot;
    AIDA::IHistogram1D *_histLeps, *_histNotLeps;
    AIDA::IHistogram1D *_histChLeps, *_histChNotLeps;
    //@}


    /// Hide the assignment operator
    MyTest& operator=(const MyTest&);
  };


}


/////////////////////////////////////


#include "Rivet/AnalysisLoader.hh"
extern "C" {
  AnalysisBuilders getAnalysisBuilders() {
    AnalysisBuilders fns;
    fns["MYTEST"] = Rivet::MyTest::create;
    return fns;
  }
}
