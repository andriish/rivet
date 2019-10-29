// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <iostream>
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class CMS_2017_I1598460 : public Analysis {
  public:

    /// Constructor
    CMS_2017_I1598460()
      : Analysis("CMS_2017_I1598460")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs,"FinalState");
      FastJets fj07(fs, FastJets::ANTIKT, 0.7);
      declare(fj07,"Jets");
      /// @todo Book histograms here, e.g.:
      for (int i=0; i<6; i++) 
	{Histo1DPtr tmp; _h_ybys.push_back(book(tmp, i+1,1,1));}
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FastJets& fj = apply<FastJets>(event,"Jets");
      const Jets& jets = fj.jetsByPt(Cuts::pt>50.*GeV && Cuts::absrap < 5.);
      // Require two jets
      if (jets.size() < 2)
        return;
      // Veto events if one of two leading jets |y|>3.0, otherwise
      // the subleading jets can become the leading jets through jet selection
      if (jets[0].absrap() > 3. || jets[1].absrap() > 3.)
        return;

      double ystar = 0.5 * std::abs(jets[0].rap() - jets[1].rap());
      double yboost = 0.5 * std::abs(jets[0].rap() + jets[1].rap());
      double ptavg = 0.5 * (jets[0].pT() + jets[1].pT());

      size_t i = (size_t)yboost;
      size_t j = (size_t)ystar;
      // index of histogram to be filled
      // yb0ys0 --> 0 yb0ys1 --> 1 ...
      size_t idx = j + 3*i - i*(i-1)/2;

      _h_ybys[idx]->fill(ptavg);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here
      for (size_t k = 0; k < _h_ybys.size(); ++k) {
        scale(_h_ybys[k], crossSection()/sumOfWeights());
      }
    }

    //@}


  private:

    //@{
    std::vector<Histo1DPtr> _h_ybys;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1598460);


}
