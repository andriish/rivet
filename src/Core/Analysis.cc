// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetYODA.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {
  Analysis::Analysis(const string& name)
    : _crossSection(-1.0),
      _gotCrossSection(false),
      _analysishandler(NULL)
  {
    ProjectionApplier::_allowProjReg = false;
    _defaultname = name;

    AnalysisInfo* ai = AnalysisInfo::make(name);
    assert(ai != 0);
    _info.reset(ai);
    assert(_info.get() != 0);
  }

  double Analysis::sqrtS() const {
    return handler().sqrtS();
  }

  const ParticlePair& Analysis::beams() const {
    return handler().beams();
  }

  const PdgIdPair Analysis::beamIds() const {
    return handler().beamIds();
  }


  const string Analysis::histoDir() const {
    /// @todo This doesn't change: calc and cache at first use!
    string path = "/" + name();
    if (handler().runName().length() > 0) {
      path = "/" + handler().runName() + path;
    }
    while (find_first(path, "//")) {
      replace_all(path, "//", "/");
    }
    return path;
  }


  const string Analysis::histoPath(const string& hname) const {
    const string path = histoDir() + "/" + hname;
    return path;
  }


  const string Analysis::histoPath(size_t datasetId, size_t xAxisId, size_t yAxisId) const {
    return histoDir() + "/" + makeAxisCode(datasetId, xAxisId, yAxisId);
  }


  const string Analysis::makeAxisCode(size_t datasetId, size_t xAxisId, size_t yAxisId) const {
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


  Log& Analysis::getLog() const {
    string logname = "Rivet.Analysis." + name();
    return Log::getLog(logname);
  }


  size_t Analysis::numEvents() const {
    return handler().numEvents();
  }


  double Analysis::sumOfWeights() const {
    return handler().sumOfWeights();
  }


  ///////////////////////////////////////////


  bool Analysis::isCompatible(const ParticlePair& beams) const {
    return isCompatible(beams.first.pdgId(),  beams.second.pdgId(),
                        beams.first.energy(), beams.second.energy());
  }


  bool Analysis::isCompatible(PdgId beam1, PdgId beam2, double e1, double e2) const {
    PdgIdPair beams(beam1, beam2);
    pair<double,double> energies(e1, e2);
    return isCompatible(beams, energies);
  }


  bool Analysis::isCompatible(const PdgIdPair& beams, const pair<double,double>& energies) const {
    // First check the beam IDs
    bool beamIdsOk = false;
    foreach (const PdgIdPair& bp, requiredBeams()) {
      if (compatible(beams, bp)) {
        beamIdsOk =  true;
        break;
      }
    }
    if (!beamIdsOk) return false;

    // Next check that the energies are compatible (within 1%, to give a bit of UI forgiveness)
    bool beamEnergiesOk = requiredEnergies().size()>0 ? false : true;
    typedef pair<double,double> DoublePair;
    foreach (const DoublePair& ep, requiredEnergies()) {
      if ((fuzzyEquals(ep.first, energies.first, 0.01) && fuzzyEquals(ep.second, energies.second, 0.01)) ||
          (fuzzyEquals(ep.first, energies.second, 0.01) && fuzzyEquals(ep.second, energies.first, 0.01))) {
        beamEnergiesOk =  true;
        break;
      }
    }
    return beamEnergiesOk;

    /// @todo Need to also check internal consistency of the analysis'
    /// beam requirements with those of the projections it uses.
  }


  ///////////////////////////////////////////


  Analysis& Analysis::setCrossSection(double xs) {
    _crossSection = xs;
    _gotCrossSection = true;
    return *this;
  }

  double Analysis::crossSection() const {
    if (!_gotCrossSection || std::isnan(_crossSection)) {
      string errMsg = "You did not set the cross section for the analysis " + name();
      throw Error(errMsg);
    }
    return _crossSection;
  }

  double Analysis::crossSectionPerEvent() const {
    const double sumW = sumOfWeights();
    assert(sumW != 0.0);
    return _crossSection / sumW;
  }



  ////////////////////////////////////////////////////////////
  // Histogramming


  void Analysis::_cacheRefData() const {
    if (_refdata.empty()) {
      MSG_TRACE("Getting refdata cache for paper " << name());
      _refdata = getRefData(name());
    }
  }


  const Scatter2D & Analysis::referenceData(const string& hname) const {
    _cacheRefData();
    MSG_TRACE("Using histo bin edges for " << name() << ":" << hname);
    return *_refdata[hname];
  }


  const Scatter2D & Analysis::referenceData(size_t datasetId, size_t xAxisId, size_t yAxisId) const {
    const string hname = makeAxisCode(datasetId, xAxisId, yAxisId);
    return referenceData(hname);
  }


  BinEdges Analysis::logBinEdges(size_t nbins, double lower, double upper) {
    assert(lower>0.0);
    assert(upper>lower);
    double loglower=log10(lower);
    double logupper=log10(upper);
    vector<double> binedges;
    double stepwidth=(logupper-loglower)/double(nbins);
    for (size_t i=0; i<=nbins; ++i) {
      binedges.push_back(pow(10.0, loglower+double(i)*stepwidth));
    }
    return binedges;
  }

  Histo1DPtr Analysis::bookHisto1D(size_t datasetId, size_t xAxisId,
				   size_t yAxisId, const string& title,
				   const string& xtitle, const string& ytitle)
  {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookHisto1D(axisCode, title, xtitle, ytitle);
  }


  Histo1DPtr Analysis::bookHisto1D(const string& hname, const string& title,
				   const string& xtitle, const string& ytitle)
  {
    // Get the bin edges (only read the AIDA file once)
    const Scatter2D & refdata = referenceData(hname);
    const string path = histoPath(hname);
    Histo1DPtr hist( new Histo1D(refdata, title) );
    addPlot(hist);
    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    // hist->setXTitle(xtitle);
    // hist->setYTitle(ytitle);
    return hist;
  }


  Histo1DPtr Analysis::bookHisto1D(const string& hname,
				   size_t nbins, double lower, double upper,
				   const string& title,
				   const string& xtitle, const string& ytitle) {
    const string path = histoPath(hname);
    Histo1DPtr hist( new Histo1D(nbins, lower, upper, path, title) );
    addPlot(hist);
    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    // hist->setXTitle(xtitle);
    // hist->setYTitle(ytitle);
    return hist;
  }


  Histo1DPtr Analysis::bookHisto1D(const string& hname,
				   const vector<double>& binedges,
				   const string& title,
				   const string& xtitle,
				   const string& ytitle) {
    const string path = histoPath(hname);
    Histo1DPtr hist( new Histo1D(binedges, path, title) );
    addPlot(hist);
    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    // hist->setXTitle(xtitle);
    // hist->setYTitle(ytitle);
    return hist;
  }

  // IHistogram2D*
  // Analysis::bookHistogram2D(const string& hname,
  // 			    size_t nxbins, double xlower, double xupper,
  // 			    size_t nybins, double ylower, double yupper,
  // 			    const string& title, const string& xtitle,
  // 			    const string& ytitle, const string& ztitle) {
  //   _makeHistoDir();
  //   const string path = histoPath(hname);
  //   IHistogram2D* hist =
  //     histogramFactory().createHistogram2D(path, title, nxbins, xlower, xupper,
  // 					   nybins, ylower, yupper);
  //   MSG_TRACE("Made 2D histogram " << hname <<  " for " << name());
  //   hist->setXTitle(xtitle);
  //   hist->setYTitle(ytitle);
  //   hist->setZTitle(ztitle);
  //   return hist;
  // }


  // IHistogram2D*
  // Analysis::bookHistogram2D(const string& hname,
  // 			    const vector<double>& xbinedges,
  // 			    const vector<double>& ybinedges,
  // 			    const string& title, const string& xtitle,
  // 			    const string& ytitle, const string& ztitle) {
  //   _makeHistoDir();
  //   const string path = histoPath(hname);
  //   IHistogram2D* hist =
  //     histogramFactory().createHistogram2D(path, title, xbinedges, ybinedges);
  //   MSG_TRACE("Made 2D histogram " << hname <<  " for " << name());
  //   hist->setXTitle(xtitle);
  //   hist->setYTitle(ytitle);
  //   hist->setZTitle(ztitle);
  //   return hist;
  // }


  /////////////////


  Profile1DPtr Analysis::bookProfile1D(size_t datasetId, size_t xAxisId,
				       size_t yAxisId, const string& title,
				       const string& xtitle, const string& ytitle) {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookProfile1D(axisCode, title, xtitle, ytitle);
  }


  Profile1DPtr Analysis::bookProfile1D(const string& hname, const string& title,
				       const string& xtitle, const string& ytitle)
  {
    // Get the bin edges (only read the AIDA file once)
    const Scatter2D & refdata = referenceData(hname);
    const string path = histoPath(hname);
    Profile1DPtr prof( new Profile1D(refdata, title) );
    addPlot(prof);
    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    // prof->setXTitle(xtitle);
    // prof->setYTitle(ytitle);
    return prof;
  }


  Profile1DPtr Analysis::bookProfile1D(const string& hname,
				       size_t nbins, double lower, double upper,
				       const string& title,
				       const string& xtitle, const string& ytitle) {
    const string path = histoPath(hname);
    Profile1DPtr prof( new Profile1D(nbins, lower, upper, path, title) );
    addPlot(prof);
    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    // prof->setXTitle(xtitle);
    // prof->setYTitle(ytitle);
    return prof;
  }


  Profile1DPtr Analysis::bookProfile1D(const string& hname,
				       const vector<double>& binedges,
				       const string& title,
				       const string& xtitle, const string& ytitle) {
    const string path = histoPath(hname);
    Profile1DPtr prof( new Profile1D(binedges, path, title) );
    addPlot(prof);
    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    // prof->setXTitle(xtitle);
    // prof->setYTitle(ytitle);
    return prof;
  }


  ///////////////////



  Scatter2DPtr Analysis::bookScatter2D(const string& hname, const string& title,
				       const string& xtitle, const string& ytitle) {
    const string path = histoPath(hname);
    Scatter2DPtr dps( new Scatter2D(path, title) );
    addPlot(dps);
    MSG_TRACE("Made data point set " << hname <<  " for " << name());
    // dps->setXTitle(xtitle);
    // dps->setYTitle(ytitle);
    return dps;
  }


  Scatter2DPtr Analysis::bookScatter2D(const string& hname,
				       size_t npts, double lower, double upper,
				       const string& title,
				       const string& xtitle, const string& ytitle) {
    Scatter2DPtr dps = bookScatter2D(hname, title, xtitle, ytitle);
    const double binwidth = (upper-lower)/npts;
    for (size_t pt = 0; pt < npts; ++pt) {
      const double bincentre = lower + (pt + 0.5) * binwidth;
      // \todo YODA check
      dps->addPoint(bincentre, 0, binwidth/2.0, 0);
      // IMeasurement* meas = dps->point(pt)->coordinate(0);
      // meas->setValue(bincentre);
      // meas->setErrorPlus(binwidth/2.0);
      // meas->setErrorMinus(binwidth/2.0);
    }
    return dps;
  }

  // \todo YODA
  // Scatter2DPtr Analysis::bookScatter2D(size_t datasetId, size_t xAxisId,
  // 				       size_t yAxisId, const string& title,
  // 				       const string& xtitle, const string& ytitle) {
  //   // Get the bin edges (only read the AIDA file once)
  //   _cacheXAxisData();
  //   // Build the axis code
  //   const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
  //   //const map<string, vector<DPSXPoint> > xpoints = getDPSXValsErrs(papername);
  //   MSG_TRACE("Using DPS x-positions for " << name() << ":" << axisCode);
  //   Scatter2DPtr dps = bookScatter2D(axisCode, title, xtitle, ytitle);
  //   const vector<Point2D> xpts = _dpsData.find(axisCode)->second;
  //   foreach ( const Point2D & pt, xpts ) {
  //     // \todo YODA check
  //     dps->addPoint(pt.x(), pt.xErrMinus(), pt.xErrPlus(), 0, 0, 0);
  //     // dps->addPoint(xpts[pt].val, xpts[pt].errminus, xpts[pt].errplus, 0, 0, 0);
  //     // IMeasurement* meas = dps->point(pt)->coordinate(0);
  //     // meas->setValue(xpts[pt].val);
  //     // meas->setErrorPlus(xpts[pt].errplus);
  //     // meas->setErrorMinus(xpts[pt].errminus);
  //   }
  //   MSG_TRACE("Made DPS " << axisCode <<  " for " << name());
  //   return dps;
  // }


  void Analysis::normalize(Histo1DPtr histo, double norm, bool includeoverflows) {
    if (!histo) {
      MSG_ERROR("Failed to normalize histo=NULL in analysis " << name() << " (norm=" << norm << ")");
      return;
    }
    MSG_TRACE("Normalizing histo " << histo->path() << " to " << norm);
    try {
      histo->normalize(norm, includeoverflows);
    } catch (YODA::WeightError& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
  }


  void Analysis::scale(Histo1DPtr histo, double scale) {
    if (!histo) {
      MSG_ERROR("Failed to scale histo=NULL in analysis " << name() << " (scale=" << scale << ")");
      return;
    }
    MSG_TRACE("Scaling histo " << histo->path() << "by factor " << scale);
    try {
      histo->scaleW(scale);
    } catch (YODA::WeightError& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
    // // Transforming the histo into a scatter after scaling
    // vector<double> x, y, ex, ey;
    // for (size_t i = 0, N = histo->numBins(); i < N; ++i) {
    //   x.push_back( histo->bin(i).midpoint() );
    //   ex.push_back(histo->bin(i).width()*0.5);
    //   y.push_back(histo->bin(i).height()*scale);
    //   ey.push_back(histo->bin(i).heightErr()*scale);
    // }
    // string title = histo->title();
    // Scatter2DPtr dps( new Scatter2D(x, y, ex, ey, hpath, title) );
    // addPlot(dps);
  }


  /// @todo 2D versions of scale and normalize... or ditch these completely?


  void Analysis::addPlot(AnalysisObjectPtr ao) {
    _plotobjects.push_back(ao);
  }

}
