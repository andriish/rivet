// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PDG_Upsilon_4S_HADRON_MULTIPLICITIES : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PDG_Upsilon_4S_HADRON_MULTIPLICITIES);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      for(int ix : _ihistos)
	_histos[ix] = bookHisto1D(ix, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double weight = event.weight();
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==300553)) {
	_wSum += weight;
	map<int,unsigned int> ncount;
	findDecayProducts(meson,ncount);
	_histos[29]->fill(10.579,double(ncount[411]+ncount[-411])*weight);
	_histos[30]->fill(10.579,double(ncount[421]+ncount[-421])*weight);
	_histos[31]->fill(10.579,double(ncount[413]+ncount[-413])*weight);
	_histos[32]->fill(10.579,double(ncount[423]+ncount[-423])*weight);
	_histos[33]->fill(10.579,double(ncount[431]+ncount[-431])*weight);
	_histos[34]->fill(10.579,double(ncount[433]+ncount[-433])*weight);
	_histos[48]->fill(10.579,double(ncount[443])*weight);
	_histos[50]->fill(10.579,double(ncount[100443])*weight);
	_histos[51]->fill(10.579,double(ncount[20443])*weight);
	_histos[53]->fill(10.579,double(ncount[445])*weight);
	_histos[60]->fill(10.579,double(ncount[321]+ncount[-321])*weight);
	_histos[61]->fill(10.579,double(ncount[321])*weight);
	_histos[62]->fill(10.579,double(ncount[-321])*weight);
	_histos[63]->fill(10.579,double(ncount[130]+ncount[310])*weight);
	_histos[64]->fill(10.579,double(ncount[323]+ncount[-323])*weight);
	_histos[65]->fill(10.579,double(ncount[313]+ncount[-313])*weight);
	_histos[87]->fill(10.579,double(ncount[211]+ncount[-211])*weight);
	_histos[88]->fill(10.579,double(ncount[111])*weight);
	_histos[89]->fill(10.579,double(ncount[221])*weight);
	_histos[90]->fill(10.579,double(ncount[113])*weight);
	_histos[92]->fill(10.579,double(ncount[333])*weight);
	_histos[96]->fill(10.579,double(ncount[4122]+ncount[-4122])*weight);
	_histos[104]->fill(10.579,double(ncount[-4222])*weight);
	_histos[106]->fill(10.579,double(ncount[-4112])*weight);
	_histos[110]->fill(10.579,double(ncount[2212]+ncount[-2212])*weight);
	_histos[113]->fill(10.579,double(ncount[3122]+ncount[-3122])*weight);
	_histos[116]->fill(10.579,double(ncount[3312]+ncount[-3312])*weight);
      }
    }
    
    void findDecayProducts(const Particle & mother,
			   map<int,unsigned int> & ncount) {
      for(const Particle & p : mother.children()) {
	int id = p.pdgId();
	if(p.children().empty()) {
	  ncount[id] += 1;
	}
	else {
	  // check particle is not a child or itself, eg copy or from photon radiation
	  bool isChild(false);
	  for(const Particle & p2 : p.children()) {
	    if(p2.pdgId()==id) {
	      isChild = true;
	      break;
	    }
	  }
	  if(!isChild) ncount[id] += 1;
	  findDecayProducts(p,ncount);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for(auto hist : _histos) {
	scale(hist.second, 1./_wSum);
      }

    }

    //@}


    /// @name Histograms
    //@{
    vector<int> _ihistos={29 ,30 ,31 ,32 ,33 ,34 ,48 ,50 ,51 ,53 ,60 ,61 ,62 ,63 ,
			  64 ,65 ,87 ,88 ,89 ,90 ,92 ,96 ,104,106,110,113,116};
    map<int,Histo1DPtr> _histos;
    double _wSum = 0.;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PDG_Upsilon_4S_HADRON_MULTIPLICITIES);


}
