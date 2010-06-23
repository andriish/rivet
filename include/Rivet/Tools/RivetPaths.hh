// -*- C++ -*-
#ifndef RIVET_RivetPaths_HH
#define RIVET_RivetPaths_HH

namespace Rivet {


  /// Get library install path
  const std::string getLibPath();

  /// Get data install path
  const std::string getDataPath();

  /// Get Rivet data install path
  const std::string getRivetDataPath();


  /// Get Rivet analysis plugin library search paths
  const std::vector<std::string> getAnalysisLibPaths();

  /// Get Rivet analysis reference data search paths
  const std::vector<std::string> getAnalysisRefPaths();

  /// Get Rivet analysis info metadata search paths
  const std::vector<std::string> getAnalysisInfoPaths();


}

#endif
