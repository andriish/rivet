#ifndef RIVET_RIVETYODA_HH
#define RIVET_RIVETYODA_HH

/// @author Andy Buckley
/// @date   2009-01-30
/// @author David Grellscheid
/// @date   2011-07-18

// Include files
#include "Rivet/Rivet.hh"
#include "Rivet/RivetYODA.fhh"
#include "YODA/AnalysisObject.h"
#include "YODA/WriterYODA.h"
#include "YODA/Histo1D.h"
#include "YODA/Profile1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Point2D.h"
#include "YODA/ReaderAIDA.h"

namespace Rivet {

  using YODA::WriterYODA;
  using YODA::ReaderAIDA;

  typedef YODA::Histo1D Histo1D;
  typedef YODA::Profile1D Profile1D;
  typedef YODA::Scatter2D Scatter2D;
  typedef YODA::Point2D Point2D;

  /// Function to get a map of all the refdata in a paper with the
  /// given @a papername.
  RefDataMap getRefData(const string& papername);

  /// Get the file system path to the AIDA reference file for this paper.
  string getDatafilePath(const string& papername);

  /// Return the integral over the histogram bins
  inline double integral(Histo1DPtr histo) {
    return histo->integral();
  }

}

#endif
