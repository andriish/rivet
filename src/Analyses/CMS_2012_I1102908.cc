// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

#include "LWH/HistogramFactory.h"
#include "LWH/Histogram1D.h"


// exclusive events are a subset of inclusive events
// so we need special error estimation treatment
void divide_set_by_subset_with_binomial_errors(AIDA::IDataPointSet* h, const AIDA::IHistogram1D & hist1, const AIDA::IHistogram1D & hist2);


#include <sstream>


namespace Rivet {


   // CMS Measurement of inclusive and exclusive dijet production
   // ratio at large rapidity intervals


class CMS_2012_I1102908 : public Analysis {
public:

  CMS_2012_I1102908() : Analysis("CMS_2012_I1102908") {}


  // ======================================================== init

  void init() {

    FinalState fs;
    FastJets akt(fs, FastJets::ANTIKT, 0.5);
    addProjection(akt, "antikT");

    // dijet ratio
    _h_dijet_ratio = bookDataPointSet(1,1,1);
    // MN dijet ratio
    _h_MN_dijet_ratio = bookDataPointSet(2,1,1);
    // exclusive dijets
    _h_DeltaY_exclusive = bookHistogram1D("DeltaY_exclusive", binEdges(1, 1, 1));
    // inclusive dijets
    _h_DeltaY_inclusive = bookHistogram1D("DeltaY_inclusive", binEdges(1, 1, 1));
    // Mueller-Navelet dijets
    _h_DeltaY_MN = bookHistogram1D("DeltaY_MN", binEdges(1, 1, 1));

  }

  // ======================================================== analyze


  void analyze(const Event & event) {
    const double weight = event.weight();

    // Jets with  pT > 35.0, -4.7 < y < 4.7
    JetAlg const &jet_alg = applyProjection<JetAlg>(event, "antikT");
    const Jets& jets = jet_alg.jets(35.0/GeV, Rivet::MAXDOUBLE, -4.7, 4.7, RAPIDITY);

    // veto event if number of jets less than 2
    if(jets.size() < 2) return;

    // loop over jet pairs
    double deltaY_MN = 0.0;
    for(Jets::const_iterator jet1 = jets.begin(); jet1 != jets.end(); ++jet1) {
      for(Jets::const_iterator jet2 = jet1 + 1; jet2 != jets.end(); ++jet2) {
        const double y1 = jet1->momentum().rapidity(), y2 = jet2->momentum().rapidity();
        const double deltaY = fabs(y1 - y2);

        // exclusive case:
        if(jets.size()==2) {
          _h_DeltaY_exclusive->fill(deltaY, weight);
        }

        // inclusive case:
        _h_DeltaY_inclusive->fill(deltaY, weight);

        // Mueller-Navelet:
        if(deltaY > deltaY_MN) {
          deltaY_MN = deltaY;
        }
      }
    }
    _h_DeltaY_MN->fill(deltaY_MN, weight);
  }

  // ======================================================== finalize


  void finalize() {
    // computing the ratio -- note that in YODA this can be replaced by a regular
    // divide with erroropt=BINOMIAL -- see comment below!!!!!
    divide_set_by_subset_with_binomial_errors(_h_dijet_ratio, *_h_DeltaY_inclusive, *_h_DeltaY_exclusive);
    divide_set_by_subset_with_binomial_errors(_h_MN_dijet_ratio, *_h_DeltaY_MN, *_h_DeltaY_exclusive);
    // removing unnecessary histograms
    histogramFactory().destroy(_h_DeltaY_inclusive);
    histogramFactory().destroy(_h_DeltaY_exclusive);
    histogramFactory().destroy(_h_DeltaY_MN);
  }

private:

  /// @name Histograms
  //@{
  AIDA::IHistogram1D*  _h_DeltaY_inclusive;
  AIDA::IHistogram1D*  _h_DeltaY_exclusive;
  AIDA::IHistogram1D*  _h_DeltaY_MN;
  AIDA::IDataPointSet* _h_dijet_ratio;
  AIDA::IDataPointSet* _h_MN_dijet_ratio;
  //@}

};

  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1102908);
}


////////// @todo YODA: This can be deleted. Use erroropt=BINOMIAL in YODA divide() instead. But we need
//////////             to be careful because this is "superset/subset" and not "subset/superset" as in the
//////////             YODA implementation, so we need to use the inverse.

// adapted from HistogramFactory::divide
// in order to change treatment of errors
void divide_set_by_subset_with_binomial_errors(AIDA::IDataPointSet* h,
                                               const AIDA::IHistogram1D & hist1, const AIDA::IHistogram1D & hist2) {
  using namespace AIDA;
  using namespace Rivet;
  using namespace LWH;
  const LWH::Histogram1D& h1 = dynamic_cast<const LWH::Histogram1D &>(hist1);
  const LWH::Histogram1D& h2 = dynamic_cast<const LWH::Histogram1D &>(hist2);
  for (int i = 0; i < h1.axis().bins(); ++i) {
    AIDA::IDataPoint* point = h->point(i);
    double yval(0), yerr(0);
    if ( h1.binHeight(i) == 0 || h2.binHeight(i) == 0 ) {
      /// @todo Bad way of handling div by zero!
      yval = 0.0;
      yerr = 0.0;
    } else {
      const double b1 = h2.binHeight(i); // subset
      const double b2 = h1.binHeight(i); // super-set
      const double w  = b1/b2;
      double e1 = h2.binError(i);
      double e2 = h1.binError(i);
      // see http://root.cern.ch/root/html/src/TH1.cxx.html#2592:
      double w_error = 0;
      if(b1 < b2) w_error = sqrt( ((1.-2.*w)*e1*e1 + w*w*e2*e2 )/(b2*b2) );
      // taking the inverse:
      yval = 1/w;
      if(w!=0) {
        yerr = w_error/pow(w,2);
      } else {
        yerr = 0;
      }
    }
    IMeasurement* y = point->coordinate(1);
    y->setValue(yval);
    y->setErrorPlus(yerr);
    y->setErrorMinus(yerr);
  }
}
