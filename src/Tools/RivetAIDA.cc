// $Id: $
// Include files 
//#include "Rivet/RivetAIDA.hh"
#include "AIDA/IAnalysisFactory.h"
#include "LWH/AnalysisFactory.h"

/// "Plugin" function to return an AIDA system (LWH impl.)
extern "C" AIDA::IAnalysisFactory* AIDA_createAnalysisFactory() {
  return new LWH::AnalysisFactory();
}
