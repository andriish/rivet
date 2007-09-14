// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
using namespace AIDA;


namespace Rivet {


  IAnalysisFactory& Analysis::analysisFactory() {
    return getHandler().analysisFactory();
  }


  ITree& Analysis::tree() {
    return getHandler().tree();
  }


  IHistogramFactory& Analysis::histogramFactory() {
    return getHandler().histogramFactory();
  }


  IDataPointSetFactory& Analysis::datapointsetFactory() {
    return getHandler().datapointsetFactory();
  }


  Log& Analysis::getLog() {
    string logname = "Rivet.Analysis." + getName();
    return Log::getLog(logname);
  }


  IHistogram1D* Analysis::bookHistogram1D(const size_t datasetId, const size_t xAxisId, 
                                          const size_t yAxisId, const string& title) {
    stringstream axisCode;
    axisCode << "ds" << datasetId << "-x" << xAxisId << "-y" << yAxisId;
    const map<string, BinEdges> data = getBinEdges(getName());
    makeHistoDir();
    const string path = getHistoDir() + "/" + axisCode.str();
    return histogramFactory().createHistogram1D(path, title, data.find(axisCode.str())->second);
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& name, const string& title, 
                                          const size_t nbins, const double lower, const double upper) {
    makeHistoDir();
    const string path = getHistoDir() + "/" + name;
    return histogramFactory().createHistogram1D(path, title, nbins, lower, upper);
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& name, const string& title, 
                                          const vector<double>& binedges) {
    makeHistoDir();
    const string path = getHistoDir() + "/" + name;
    return histogramFactory().createHistogram1D(path, title, binedges);
  }


  IDataPointSet* Analysis::bookDataPointSet(const string& name, const string& title) {
    makeHistoDir();
    const string path = getHistoDir() + "/" + name;
    return datapointsetFactory().create(path, title, 2);
  }


  IDataPointSet* Analysis::bookDataPointSet(const string& name, const string& title, 
                                            const size_t npts, const double lower, const double upper) {
    IDataPointSet* dps = bookDataPointSet(name, title);
    for (size_t pt = 0; pt < npts; ++pt) {
      const double binwidth = (upper-lower)/npts;
      const double bincentre = lower + (pt + 0.5) * binwidth;
      dps->addPoint();
      IMeasurement* meas = dps->point(pt)->coordinate(0);
      meas->setValue(bincentre);
      meas->setErrorPlus(binwidth/2.0);
      meas->setErrorMinus(binwidth/2.0);
    }
    return dps;
  }


  IDataPointSet* Analysis::bookDataPointSet(const size_t datasetId, const size_t xAxisId, 
                                            const size_t yAxisId, const string& title) {
    /// @todo Implement this?
    throw runtime_error("Auto-booking of DataPointSets is not yet implemented");
  }


  inline void Analysis::makeHistoDir() {
    if (!_madeHistoDir) {
      if (! getName().empty()) tree().mkdir(getHistoDir());
      _madeHistoDir = true;
    }
  }


  const Cuts Analysis::getCuts() const {
    Cuts totalCuts = _cuts;
    for (set<Projection*>::const_iterator p = _projections.begin(); p != _projections.end(); ++p) {
      totalCuts.addCuts((*p)->getCuts());
    }
    return totalCuts;
  }
  

  const bool Analysis::checkConsistency() const {
    // Check consistency of analysis beams with allowed beams of each contained projection.
    for (set<Projection*>::const_iterator p = _projections.begin(); p != _projections.end(); ++p) {
      if (! compatible(getBeams(), (*p)->getBeamPairs()) ) {
        throw runtime_error("Analysis " + getName() + " beams are inconsistent with " 
                            + "allowed beams for projection " + (*p)->getName());
      }
    }
    // Check the consistency of the accumulated cuts (throws if wrong).
    getCuts().checkConsistency();
    return true;
  }


  set<Projection*> Analysis::getProjections() const {
    set<Projection*> totalProjections = _projections;
    for (set<Projection*>::const_iterator p = _projections.begin(); p != _projections.end(); ++p) {
      totalProjections.insert((*p)->getProjections().begin(), (*p)->getProjections().end());
    }
    return totalProjections;
  }


}
