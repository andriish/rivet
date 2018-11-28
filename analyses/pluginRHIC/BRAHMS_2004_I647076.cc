// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  /// @brief BRAHMS Centrality projection.
  class BRAHMSCentrality : public SingleValueProjection {
  public:
    // Constructor
    BRAHMSCentrality() : SingleValueProjection() {
      setName("BRAHMSCentrality");
      // Using here the BRAHMS reaction centrality from eg. 1602.01183, which
      // might not be correct.
      const ChargedFinalState cfs(Cuts::pT > 0.1*GeV && Cuts::abseta < 2.2);
      this->declare(cfs,"ChargedFinalState");
    }
    // Destructor
    virtual ~BRAHMSCentrality() {}

    // Do the projection. Count the number of charged particles in
    // the specified range.
    virtual void project(const Event& e) {
      clear();
      set(apply<ChargedFinalState>
        (e, "ChargedFinalState").particles().size());
    }

    // Compare to another projection.
    virtual int compare(const Projection& p) const {
      // This projection is only used for the analysis below.
      return UNDEFINED;
    }

    // Clone this projection
    virtual unique_ptr<Projection> clone() const {
      return unique_ptr<Projection>(new BRAHMSCentrality(*this));
    }
  };

  /// @brief  Brahms centrality analysis
  class BRAHMS_2004_CENTRALITY : public Analysis {
  public:
    // Constructor
    BRAHMS_2004_CENTRALITY() : Analysis("BRAHMS_2004_CENTRALITY") {}

    // Initialize the analysis
    void init() {
       declare(BRAHMSCentrality(),"BCEN");
       mult = bookHisto1D("mult",450,0,4500);
       imp = bookHisto1D("imp",100,0,20);
    }

    // Analyse a single event
    void analyze(const Event& event) {
      // Get and fill in the impact parameter value if the information
      // is valid.
      const HepMC::GenEvent* ge = event.genEvent();
      const HepMC::HeavyIon* hi = ge->heavy_ion();
      if (hi && hi->is_valid())
	imp->fill(hi->impact_parameter(), event.weight());
    
      mult->fill(apply<BRAHMSCentrality>(event,"BCEN")(), event.weight());
    }

    // Finalize the analysis
    void finalize() {
      if(mult->numEntries()) mult->normalize();
      if(imp->numEntries()) imp->normalize();
    
    }

  private:
    // Histograms.
    Histo1DPtr mult;
    Histo1DPtr imp;
  
  };
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BRAHMS_2004_CENTRALITY);

  /// @brief Brahms pT spectra for id particles AuAu @ 200GeV/nn
  class BRAHMS_2004_I647076 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BRAHMS_2004_I647076);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // Centrality Projection
      declareCentrality(BRAHMSCentrality(), "BRAHMS_2004_Centrality","mult","mult");
      // TODO: Change to PrimaryParticles depending on feed down.
      declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");

      rapIntervals = {{-0.1,0.},{0.,0.1}};
      // Book histograms
      for (int i = 1, N = rapIntervals.size(); i <= N; ++i) {
        piPlus.push_back(bookHisto1D(1, 1, i));
      }
      centSow = bookCounter("centSow");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double w = event.weight();
      const CentralityProjection& cent = apply<CentralityProjection>(event,"mult");
      if (cent() > 5.0) return;
      centSow->fill(w);
      const FinalState& fs = apply<FinalState>(event,"FS");
      for (const auto& p : fs.particles()) {
        double y = p.rapidity();
	for (int i = 0, N = rapIntervals.size(); i < N; ++i) {
	  if (y > rapIntervals[i].first && y <= rapIntervals[i].second) {
	    double dy = rapIntervals[i].second - rapIntervals[i].first;
	    double pT = p.pT();
	    double nWeight = w / ( 2*M_PI*pT*dy);
	    if (p.pid() == 211) piPlus[i]->fill(pT, nWeight);
	  }
	    
	} 
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (int i = 0, N = rapIntervals.size(); i < N; ++i) {
        piPlus[i]->scaleW(1./centSow->sumW());
      }

    }

    //@}

    // The rapidity intervals.
    vector<pair<double, double> > rapIntervals;

    /// @name Histograms
    //@{
    vector<Histo1DPtr> piPlus;
    CounterPtr centSow;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BRAHMS_2004_I647076);


}
