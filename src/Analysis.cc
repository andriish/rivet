// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "AIDA/IManagedObject.h"
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


  size_t Analysis::numEvents() { return getHandler().numEvents(); }


  double Analysis::sumOfWeights() { return getHandler().sumOfWeights(); }


  void Analysis::_cacheBinEdges() {
    if (_histBinEdges.empty()) {
      getLog() << Log::TRACE << "Getting histo bin edges from AIDA for paper " << getName() << endl;
      _histBinEdges = getBinEdges(getName());
    }
  }


  string Analysis::_makeAxisCode(const size_t datasetId, const size_t xAxisId, const size_t yAxisId) {
    stringstream axisCode;
    axisCode << "d";
    if (datasetId < 10) axisCode << 0;
    axisCode << datasetId;
    axisCode << "-x";
    if (xAxisId < 10) axisCode << 0;
    axisCode << xAxisId;
    axisCode << "-y";
    if (yAxisId < 10) axisCode << 0;
    axisCode << yAxisId;
    return axisCode.str();
  }


  IHistogram1D* Analysis::bookHistogram1D(const size_t datasetId, const size_t xAxisId, 
                                          const size_t yAxisId, const string& title) {
    // Get the bin edges (only read the AIDA file once)
    _cacheBinEdges();
    // Build the axis code
    const string axisCode = _makeAxisCode(datasetId, xAxisId, yAxisId);
    getLog() << Log::TRACE << "Using histo bin edges for " << getName() << ":" << axisCode << endl;
    const BinEdges edges = _histBinEdges.find(axisCode)->second;
    _makeHistoDir();
    const string path = getHistoDir() + "/" + axisCode;
    IHistogram1D* hist = histogramFactory().createHistogram1D(path, title, edges);
    getLog() << Log::TRACE << "Made histogram " << axisCode <<  " for " << getName() << endl;
    return hist;
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& name, const string& title, 
                                          const size_t nbins, const double lower, const double upper) {
    _makeHistoDir();
    const string path = getHistoDir() + "/" + name;
    return histogramFactory().createHistogram1D(path, title, nbins, lower, upper);
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& name, const string& title, 
                                          const vector<double>& binedges) {
    _makeHistoDir();
    const string path = getHistoDir() + "/" + name;
    return histogramFactory().createHistogram1D(path, title, binedges);
  }



  IProfile1D* Analysis::bookProfile1D(const size_t datasetId, const size_t xAxisId, 
                                          const size_t yAxisId, const string& title) {
    // Get the bin edges (only read the AIDA file once)
    _cacheBinEdges();
    // Build the axis code
    const string axisCode = _makeAxisCode(datasetId, xAxisId, yAxisId);
    getLog() << Log::TRACE << "Using profile histo bin edges for " << getName() << ":" << axisCode << endl;
    const BinEdges edges = _histBinEdges.find(axisCode)->second;
    _makeHistoDir();
    const string path = getHistoDir() + "/" + axisCode;
    IProfile1D* prof = histogramFactory().createProfile1D(path, title, edges);
    getLog() << Log::TRACE << "Made profile histogram " << axisCode <<  " for " << getName() << endl;
    return prof;
  }


  IProfile1D* Analysis::bookProfile1D(const string& name, const string& title, 
                                          const size_t nbins, const double lower, const double upper) {
    _makeHistoDir();
    const string path = getHistoDir() + "/" + name;
    return histogramFactory().createProfile1D(path, title, nbins, lower, upper);
  }


  IProfile1D* Analysis::bookProfile1D(const string& name, const string& title, 
                                          const vector<double>& binedges) {
    _makeHistoDir();
    const string path = getHistoDir() + "/" + name;
    return histogramFactory().createProfile1D(path, title, binedges);
  }


  IDataPointSet* Analysis::bookDataPointSet(const string& name, const string& title) {
    _makeHistoDir();
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


  inline void Analysis::_makeHistoDir() {
    if (!_madeHistoDir) {
      if (! getName().empty()) tree().mkdir(getHistoDir());
      _madeHistoDir = true;
    }
  }


  const Cuts Analysis::getCuts() const {
    Cuts totalCuts = _cuts;
    for (set<ConstProjectionPtr>::const_iterator p = _projections.begin(); p != _projections.end(); ++p) {
      totalCuts.addCuts((*p)->getCuts());
    }
    return totalCuts;
  }
  

  const bool Analysis::checkConsistency() const {
    // Check consistency of analysis beams with allowed beams of each contained projection.
    for (set<ConstProjectionPtr>::const_iterator p = _projections.begin(); p != _projections.end(); ++p) {
      if (! compatible(getBeams(), (*p)->getBeamPairs()) ) {
        throw runtime_error("Analysis " + getName() + " beams are inconsistent with " 
                            + "allowed beams for projection " + (*p)->getName());
      }
    }
    // Check the consistency of the accumulated cuts (throws if wrong).
    getCuts().checkConsistency();
    return true;
  }


  set<ConstProjectionPtr> Analysis::getProjections() const {
    set<ConstProjectionPtr> totalProjections = _projections;
    for (set<ConstProjectionPtr>::const_iterator p = _projections.begin(); p != _projections.end(); ++p) {
      totalProjections.insert((*p)->getProjections().begin(), (*p)->getProjections().end());
    }
    return totalProjections;
  }

void Analysis::normalize(AIDA::IHistogram1D*& histo, const double norm) {
  // Calculate histogram area even with non-uniform bin widths
  // (sumAllBinHeights works for uniform bin widths)
  Log& log = getLog();
  log << Log::TRACE << "Normalizing histo " << histo->title() << " to " << norm << endl;
  
  double oldintg = 0.0;
  int nBins = histo->axis().bins();
  for(int iBin = 0; iBin != nBins; ++iBin){
    /// @todo Leaving out factor of binWidth because AIDA's "height" already includes a width factor
    oldintg += histo->binHeight(iBin);// * histo->axis().binWidth(iBin);
  }
  if (oldintg == 0.0) return;

  const double scale = norm/oldintg;
  std::vector<double> x, y, ex, ey;
  for ( int i = 0, N = histo->axis().bins(); i < N; ++i ) {
    x.push_back(0.5 * (histo->axis().binLowerEdge(i) + histo->axis().binUpperEdge(i)));
    ex.push_back(histo->axis().binWidth(i)*0.5);
    ///@todo "Bin height" is a misnomer in the AIDA spec: width is neglected
    y.push_back(histo->binHeight(i)*scale/histo->axis().binWidth(i));
    // We'd like to do this: y.push_back(histo->binHeight(i) * scale);
    ///@todo "Bin error" is a misnomer in the AIDA spec: width is neglected
    ey.push_back(histo->binError(i)*scale/(0.5*histo->axis().binWidth(i)));
    // We'd like to do this: ey.push_back(histo->binError(i) * scale);
  }

  std::string path =
    tree().findPath(dynamic_cast<AIDA::IManagedObject&>(*histo));
  std::string title = histo->title();
  tree().mkdir("/tmpnormalize");
  tree().mv(path, "/tmpnormalize");

  datapointsetFactory().createXY(path, title, x, y, ex, ey);

  tree().rm(tree().findPath(dynamic_cast<AIDA::IManagedObject&>(*histo)));
  tree().rmdir("/tmpnormalize");
  histo = 0;

}


}
