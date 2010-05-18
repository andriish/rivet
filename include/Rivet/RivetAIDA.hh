#ifndef RIVET_RIVETAIDA_HH
#define RIVET_RIVETAIDA_HH

/// @author Andy Buckley
/// @date   2009-01-30

// Include files
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.fhh"
#include "LWH/AIAnalysisFactory.h"
#include "LWH/AIHistogramFactory.h"
#include "LWH/AIHistogram1D.h"
#include "LWH/AIProfile1D.h"
#include "LWH/AITreeFactory.h"
#include "LWH/AIDataPointSetFactory.h"
#include "LWH/AIDataPointSet.h"
#include "LWH/AIDataPoint.h"
#include "LWH/AIMeasurement.h"
#include "LWH/AITree.h"
#include "LWH/AIAxis.h"


namespace Rivet {

  AIDA::IAnalysisFactory* createAnalysisFactory();

  /// Function to get a map of all the bin edge vectors in a paper with the
  /// given @a papername.
  const map<string, BinEdges> getBinEdges(string papername);

  const map<string, BinEdges>
  getBinEdges(const map<string, vector<DPSXPoint> >& xpoints);

  const map<string, vector<DPSXPoint> > getDPSXValsErrs(string papername);

  /// Get the file system path to the AIDA reference file for this paper.
  const string getDataPath(string papername);

  /// Return the integral over the histogram bins assuming it has been
  // normalize()d.
  inline double integral(AIDA::IHistogram1D* histo) {
    double intg = 0.;
    for ( int i = 0; i < histo->axis().bins(); ++i )
      // Don't multiply with binWidth -- it's already included in binHeight
      intg += histo->binHeight(i); // * histo->axis().binWidth(i);
    return intg;
  }



  using AIDA::IHistogram1D;
  using AIDA::IDataPointSet;
  using AIDA::IDataPoint;
  using AIDA::IMeasurement;
  using AIDA::ITree;
  using AIDA::IAxis;
  using AIDA::IProfile1D;

}

#endif
