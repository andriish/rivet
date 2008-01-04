#ifndef RIVETAIDA_HH 
#define RIVETAIDA_HH 

/// @author Andy Buckley
/// @date   2007-11-08

// Include files
#include "Rivet/Rivet.hh"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IProfile1D.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/IDataPointSetFactory.h"
#include "AIDA/IDataPointSet.h"
#include "AIDA/IDataPoint.h"
#include "AIDA/IMeasurement.h"
#include "AIDA/ITree.h"
#include "AIDA/IAxis.h"


namespace Rivet {

  /// Function to get a map of all the bin edge vectors in a paper with the
  /// given @a papername.
  const map<string, BinEdges> getBinEdges(string papername);

  /// Get the file system path to the AIDA reference file for this paper.
  const string getDataPath(string papername);

  /// Normalize the histogram to @a norm .
  inline void normalize(AIDA::IHistogram1D* histo, const double norm=1.0) {
    assert(norm != 0.0);
    double area = 0;
    for (int i=0; i < histo->axis().bins(); ++i) {
      area += histo->binHeight(i) * histo->axis().binWidth(i);
    }
    if (area != 0) {
      histo->scale(norm/area);
    }
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
