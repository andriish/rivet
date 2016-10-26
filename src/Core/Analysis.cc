// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/BeamConstraint.hh"

namespace Rivet {


  Analysis::Analysis(const string& name)
    : _crossSection(-1.0),
      _gotCrossSection(false),
      _analysishandler(NULL)
  {
    ProjectionApplier::_allowProjReg = false;
    _defaultname = name;

    unique_ptr<AnalysisInfo> ai = AnalysisInfo::make(name);
    assert(ai);
    _info = move(ai);
    assert(_info);
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
    /// @todo Cache in a member variable
    string _histoDir;
    if (_histoDir.empty()) {
      _histoDir = "/" + name();
      if (handler().runName().length() > 0) {
        _histoDir = "/" + handler().runName() + _histoDir;
      }
      replace_all(_histoDir, "//", "/"); //< iterates until none
    }
    return _histoDir;
  }


  const string Analysis::histoPath(const string& hname) const {
    const string path = histoDir() + "/" + hname;
    return path;
  }


  const string Analysis::histoPath(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    return histoDir() + "/" + makeAxisCode(datasetId, xAxisId, yAxisId);
  }


  const string Analysis::makeAxisCode(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
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
    return isCompatible(beams.first.pid(),  beams.second.pid(),
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

    // Next check that the energies are compatible (within 1% or 1 GeV, whichever is larger, for a bit of UI forgiveness)

    /// @todo Use some sort of standard ordering to improve comparisons, esp. when the two beams are different particles
    bool beamEnergiesOk = requiredEnergies().size() > 0 ? false : true;
    typedef pair<double,double> DoublePair;
    foreach (const DoublePair& ep, requiredEnergies()) {
      if ((fuzzyEquals(ep.first, energies.first, 0.01) && fuzzyEquals(ep.second, energies.second, 0.01)) ||
          (fuzzyEquals(ep.first, energies.second, 0.01) && fuzzyEquals(ep.second, energies.first, 0.01)) ||
          (abs(ep.first - energies.first) < 1*GeV && abs(ep.second - energies.second) < 1*GeV) ||
          (abs(ep.first - energies.second) < 1*GeV && abs(ep.second - energies.first) < 1*GeV)) {
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
    return _crossSection/sumOfWeights();
  }



  ////////////////////////////////////////////////////////////
  // Histogramming


  void Analysis::_cacheRefData() const {
    if (_refdata.empty()) {
      MSG_TRACE("Getting refdata cache for paper " << name());
      _refdata = getRefData(name());
    }
  }


  const Scatter2D& Analysis::refData(const string& hname) const {
    _cacheRefData();
    MSG_TRACE("Using histo bin edges for " << name() << ":" << hname);
    if (!_refdata[hname]) {
      MSG_ERROR("Can't find reference histogram " << hname);
      throw Exception("Reference data " + hname + " not found.");
    }
    return dynamic_cast<Scatter2D&>(*_refdata[hname]);
  }


  const Scatter2D& Analysis::refData(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    const string hname = makeAxisCode(datasetId, xAxisId, yAxisId);
    return refData(hname);
  }


  CounterPtr& Analysis::bookCounter(const string& cname,
                                   const string& title) {
    const string path = histoPath(cname);
    shared_ptr<CounterPtr> ctr =
        make_shared<CounterPtr>(handler().numWeights(), Counter(path, title));
    addAnalysisObject(ctr);
    MSG_TRACE("Made counter " << cname << " for " << name());
    return *ctr;
  }


  CounterPtr& Analysis::bookCounter(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                   const string& title) {
                                   // const string& xtitle,
                                   // const string& ytitle) {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookCounter(axisCode, title);
  }


  Histo1DPtr& Analysis::bookHisto1D(const string& hname,
                                    size_t nbins, double lower, double upper,
                                    const string& title,
                                    const string& xtitle,
                                    const string& ytitle) {
    const string path = histoPath(hname);

    Histo1D hist = Histo1D(nbins, lower, upper, path, title);

    hist.setAnnotation("XLabel", xtitle);
    hist.setAnnotation("YLabel", ytitle);

    shared_ptr<Histo1DPtr> hp =
        make_shared<Histo1DPtr>(handler().numWeights(), hist);

    addAnalysisObject(hp);
    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    return *hp;
  }


  Histo1DPtr& Analysis::bookHisto1D(const string& hname,
                                   const vector<double>& binedges,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    const string path = histoPath(hname);

    Histo1D hist = Histo1D(binedges, path, title);

    hist.setAnnotation("XLabel", xtitle);
    hist.setAnnotation("YLabel", ytitle);


    shared_ptr<Histo1DPtr> hp =
        make_shared<Histo1DPtr>(handler().numWeights(), hist);
    addAnalysisObject(hp);

    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    return *hp;
  }


  Histo1DPtr& Analysis::bookHisto1D(const string& hname,
                                   const Scatter2D& refscatter,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    const string path = histoPath(hname);
    Histo1D hist = Histo1D(refscatter, path);

    hist.setTitle(title);
    hist.setAnnotation("XLabel", xtitle);
    hist.setAnnotation("YLabel", ytitle);

    shared_ptr<Histo1DPtr> hp = 
        make_shared<Histo1DPtr>(handler().numWeights(), hist);
    addAnalysisObject(hp);

    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    return *hp;
  }


  Histo1DPtr& Analysis::bookHisto1D(const string& hname,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    const Scatter2D& refdata = refData(hname);
    return bookHisto1D(hname, refdata, title, xtitle, ytitle);
  }


  Histo1DPtr& Analysis::bookHisto1D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookHisto1D(axisCode, title, xtitle, ytitle);
  }


  /// @todo Add booking methods which take a path, titles and *a reference Scatter from which to book*


  /////////////////


  Histo2DPtr& Analysis::bookHisto2D(const string& hname,
                                   size_t nxbins, double xlower, double xupper,
                                   size_t nybins, double ylower, double yupper,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle)
  {
    const string path = histoPath(hname);

    Histo2D hist(nxbins, xlower, xupper, nybins, ylower, yupper, path, title);
    hist.setAnnotation("XLabel", xtitle);
    hist.setAnnotation("YLabel", ytitle);
    hist.setAnnotation("ZLabel", ztitle);

    shared_ptr<Histo2DPtr> hp =
        make_shared<Histo2DPtr>(handler().numWeights(), hist);
    addAnalysisObject(hp);

    MSG_TRACE("Made 2D histogram " << hname <<  " for " << name());
    return *hp;
  }


  Histo2DPtr& Analysis::bookHisto2D(const string& hname,
                                   const vector<double>& xbinedges,
                                   const vector<double>& ybinedges,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle)
  {
    const string path = histoPath(hname);

    Histo2D hist(xbinedges, ybinedges, path, title);
    hist.setAnnotation("XLabel", xtitle);
    hist.setAnnotation("YLabel", ytitle);
    hist.setAnnotation("ZLabel", ztitle);

    shared_ptr<Histo2DPtr> hp =
        make_shared<Histo2DPtr>(handler().numWeights(), hist);
    addAnalysisObject(hp);

    MSG_TRACE("Made 2D histogram " << hname <<  " for " << name());
    return *hp;
  }


  Profile1DPtr& Analysis::bookProfile1D(const string& hname,
                                       size_t nbins, double lower, double upper,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);

    Profile1D prof(nbins, lower, upper, path, title);
    prof.setAnnotation("XLabel", xtitle);
    prof.setAnnotation("YLabel", ytitle);

    shared_ptr<Profile1DPtr> pp =
        make_shared<Profile1DPtr>(handler().numWeights(), prof);
    addAnalysisObject(pp);

    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    return *pp;
  }


  Profile1DPtr& Analysis::bookProfile1D(const string& hname,
                                       const vector<double>& binedges,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);
    Profile1D prof(binedges, path, title);
    prof.setAnnotation("XLabel", xtitle);
    prof.setAnnotation("YLabel", ytitle);

    shared_ptr<Profile1DPtr> pp =
        make_shared<Profile1DPtr>(handler().numWeights(), prof);
    addAnalysisObject(pp);

    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    return *pp;
  }


  Profile1DPtr& Analysis::bookProfile1D(const string& hname,
                                       const Scatter2D& refscatter,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);
    Profile1D prof(refscatter, path);
    prof.setTitle(title);
    prof.setAnnotation("XLabel", xtitle);
    prof.setAnnotation("YLabel", ytitle);

    shared_ptr<Profile1DPtr> pp =
        make_shared<Profile1DPtr>(handler().numWeights(), prof);
    addAnalysisObject(pp);

    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    return *pp;
  }


  Profile1DPtr& Analysis::bookProfile1D(const string& hname,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const Scatter2D& refdata = refData(hname);
    return bookProfile1D(hname, refdata, title, xtitle, ytitle);
  }


  Profile1DPtr& Analysis::bookProfile1D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookProfile1D(axisCode, title, xtitle, ytitle);
  }


  Profile2DPtr& Analysis::bookProfile2D(const string& hname,
                                   size_t nxbins, double xlower, double xupper,
                                   size_t nybins, double ylower, double yupper,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle)
  {
    const string path = histoPath(hname);
    Profile2D prof(nxbins, xlower, xupper, nybins, ylower, yupper, path, title);
    prof.setAnnotation("XLabel", xtitle);
    prof.setAnnotation("YLabel", ytitle);
    prof.setAnnotation("ZLabel", ztitle);

    shared_ptr<Profile2DPtr> pp =
        make_shared<Profile2DPtr>(handler().numWeights(), prof);
    addAnalysisObject(pp);

    MSG_TRACE("Made 2D profile histogram " << hname <<  " for " << name());
    return *pp;
  }


  Profile2DPtr& Analysis::bookProfile2D(const string& hname,
                                   const vector<double>& xbinedges,
                                   const vector<double>& ybinedges,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle)
  {
    const string path = histoPath(hname);
    Profile2D prof(xbinedges, ybinedges, path, title);
    prof.setAnnotation("XLabel", xtitle);
    prof.setAnnotation("YLabel", ytitle);
    prof.setAnnotation("ZLabel", ztitle);

    shared_ptr<Profile2DPtr> pp =
        make_shared<Profile2DPtr>(handler().numWeights(), prof);
    addAnalysisObject(pp);

    MSG_TRACE("Made 2D profile histogram " << hname <<  " for " << name());
    return *pp;
  }


  Scatter2DPtr& Analysis::bookScatter2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                       bool copy_pts,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookScatter2D(axisCode, copy_pts, title, xtitle, ytitle);
  }


  Scatter2DPtr& Analysis::bookScatter2D(const string& hname,
                                       bool copy_pts,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    Scatter2D scat;
    const string path = histoPath(hname);
    if (copy_pts) {
      const Scatter2D& refdata = refData(hname);
      scat = Scatter2D(refdata, path);
      foreach (Point2D& p, scat.points()) p.setY(0, 0);
    } else {
      scat = Scatter2D(path);
    }

    scat.setTitle(title);
    scat.setAnnotation("XLabel", xtitle);
    scat.setAnnotation("YLabel", ytitle);

    shared_ptr<Scatter2DPtr> sp = make_shared<Scatter2DPtr>(scat);
    addAnalysisObject(sp);

    MSG_TRACE("Made scatter " << hname <<  " for " << name());
    return *sp;
  }


  Scatter2DPtr& Analysis::bookScatter2D(const string& hname,
                                       size_t npts, double lower, double upper,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);
    Scatter2D scat;
    const double binwidth = (upper-lower)/npts;
    for (size_t pt = 0; pt < npts; ++pt) {
      const double bincentre = lower + (pt + 0.5) * binwidth;
      scat.addPoint(bincentre, 0, binwidth/2.0, 0);
    }
    scat.setTitle(title);
    scat.setAnnotation("XLabel", xtitle);
    scat.setAnnotation("YLabel", ytitle);

    shared_ptr<Scatter2DPtr> sp = make_shared<Scatter2DPtr>(scat);
    addAnalysisObject(sp);

    MSG_TRACE("Made scatter " << hname <<  " for " << name());
    return *sp;
  }


  Scatter2DPtr& Analysis::bookScatter2D(const string& hname,
                                       const vector<double>& binedges,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);
    Scatter2D scat;
    for (size_t pt = 0; pt < binedges.size()-1; ++pt) {
      const double bincentre = (binedges[pt] + binedges[pt+1]) / 2.0;
      const double binwidth = binedges[pt+1] - binedges[pt];
      scat.addPoint(bincentre, 0, binwidth/2.0, 0);
    }

    scat.setTitle(title);
    scat.setAnnotation("XLabel", xtitle);
    scat.setAnnotation("YLabel", ytitle);

    shared_ptr<Scatter2DPtr> sp = make_shared<Scatter2DPtr>(scat);

    MSG_TRACE("Made scatter " << hname <<  " for " << name());
    return *sp;
  }


  void Analysis::divide(CounterPtr c1, CounterPtr c2, Scatter1DPtr s) const {
    const string path = s->path();
    *s = *c1 / *c2;
    s->setPath(path);
  }

  void Analysis::divide(const Counter& c1, const Counter& c2, Scatter1DPtr s) const {
    const string path = s->path();
    *s = c1 / c2;
    s->setPath(path);
  }


  void Analysis::divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = *h1 / *h2;
    s->setPath(path);
  }

  void Analysis::divide(const Histo1D& h1, const Histo1D& h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = h1 / h2;
    s->setPath(path);
  }


  void Analysis::divide(Profile1DPtr p1, Profile1DPtr p2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = *p1 / *p2;
    s->setPath(path);
  }

  void Analysis::divide(const Profile1D& p1, const Profile1D& p2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = p1 / p2;
    s->setPath(path);
  }


  void Analysis::divide(Histo2DPtr h1, Histo2DPtr h2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = *h1 / *h2;
    s->setPath(path);
  }

  void Analysis::divide(const Histo2D& h1, const Histo2D& h2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = h1 / h2;
    s->setPath(path);
  }


  void Analysis::divide(Profile2DPtr p1, Profile2DPtr p2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = *p1 / *p2;
    s->setPath(path);
  }

  void Analysis::divide(const Profile2D& p1, const Profile2D& p2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = p1 / p2;
    s->setPath(path);
  }


  /// @todo Counter and Histo2D efficiencies and asymms


  void Analysis::efficiency(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::efficiency(*h1, *h2);
    s->setPath(path);
  }

  void Analysis::efficiency(const Histo1D& h1, const Histo1D& h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::efficiency(h1, h2);
    s->setPath(path);
  }


  void Analysis::asymm(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::asymm(*h1, *h2);
    s->setPath(path);
  }

  void Analysis::asymm(const Histo1D& h1, const Histo1D& h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::asymm(h1, h2);
    s->setPath(path);
  }


  void Analysis::scale(CounterPtr cnt, double factor) {
    if (!cnt) {
      MSG_WARNING("Failed to scale counter=NULL in analysis " << name() << " (scale=" << factor << ")");
      return;
    }
    if (std::isnan(factor) || std::isinf(factor)) {
      MSG_WARNING("Failed to scale counter=" << cnt->path() << " in analysis: " << name() << " (invalid scale factor = " << factor << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling counter " << cnt->path() << " by factor " << factor);
    try {
      cnt->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale counter " << cnt->path());
      return;
    }
  }


  void Analysis::normalize(Histo1DPtr histo, double norm, bool includeoverflows) {
    if (!histo) {
      MSG_WARNING("Failed to normalize histo=NULL in analysis " << name() << " (norm=" << norm << ")");
      return;
    }
    MSG_TRACE("Normalizing histo " << histo->path() << " to " << norm);
    try {
      histo->normalize(norm, includeoverflows);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
  }


  void Analysis::scale(Histo1DPtr histo, double factor) {
    if (!histo) {
      MSG_WARNING("Failed to scale histo=NULL in analysis " << name() << " (scale=" << factor << ")");
      return;
    }
    if (std::isnan(factor) || std::isinf(factor)) {
      MSG_WARNING("Failed to scale histo=" << histo->path() << " in analysis: " << name() << " (invalid scale factor = " << factor << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling histo " << histo->path() << " by factor " << factor);
    try {
      histo->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale histo " << histo->path());
      return;
    }
  }


  void Analysis::normalize(Histo2DPtr histo, double norm, bool includeoverflows) {
    if (!histo) {
      MSG_ERROR("Failed to normalize histo=NULL in analysis " << name() << " (norm=" << norm << ")");
      return;
    }
    MSG_TRACE("Normalizing histo " << histo->path() << " to " << norm);
    try {
      histo->normalize(norm, includeoverflows);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
  }


  void Analysis::scale(Histo2DPtr histo, double factor) {
    if (!histo) {
      MSG_ERROR("Failed to scale histo=NULL in analysis " << name() << " (scale=" << factor << ")");
      return;
    }
    if (std::isnan(factor) || std::isinf(factor)) {
      MSG_ERROR("Failed to scale histo=" << histo->path() << " in analysis: " << name() << " (invalid scale factor = " << factor << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling histo " << histo->path() << " by factor " << factor);
    try {
      histo->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale histo " << histo->path());
      return;
    }
  }


  void Analysis::integrate(Histo1DPtr h, Scatter2DPtr s) const {
    // preserve the path info
    const string path = s->path();
    *s = toIntegralHisto(*h);
    s->setPath(path);
  }

  void Analysis::integrate(const Histo1D& h, Scatter2DPtr s) const {
    // preserve the path info
    const string path = s->path();
    *s = toIntegralHisto(h);
    s->setPath(path);
  }


  /// @todo 2D versions of integrate... defined how, exactly?!?


  //////////////////////////////////


  // @todo
  // special handling for scatters
  void Analysis::addAnalysisObject(const shared_ptr<Scatter1DPtr>& ao) {
    _scatters.push_back(ao);
  }
  void Analysis::addAnalysisObject(const shared_ptr<Scatter2DPtr>& ao) {
    _scatters.push_back(ao);
  }
  void Analysis::addAnalysisObject(const shared_ptr<Scatter3DPtr>& ao) {
    _scatters.push_back(ao);
  }

  void Analysis::addAnalysisObject(const shared_ptr<MultiweightAOPtr>& ao) {
    _analysisobjects.push_back(ao);
  }

  void Analysis::removeAnalysisObject(const string& path) {
    for (vector<shared_ptr<MultiweightAOPtr> >::iterator it = _analysisobjects.begin();  it != _analysisobjects.end(); ++it) {
      if ((**it)->path() == path) {
        _analysisobjects.erase(it);
        break;
      }
    }

    for (vector<shared_ptr<AnalysisObjectPtr> >::iterator it = _scatters.begin();  it != _scatters.end(); ++it) {
      if ((**it)->path() == path) {
        _scatters.erase(it);
        break;
      }
    }
  }

  /// @todo can we really remove (multiweighted) analysis objects by == operator??
  void Analysis::removeAnalysisObject(const Scatter1DPtr& ao) {
    for (vector<shared_ptr<AnalysisObjectPtr> >::iterator it = _scatters.begin();  it != _scatters.end(); ++it) {
      if (**it == ao) {
        _scatters.erase(it);
        break;
      }
    }
  }

  void Analysis::removeAnalysisObject(const Scatter2DPtr& ao) {
    for (vector<shared_ptr<AnalysisObjectPtr> >::iterator it = _scatters.begin();  it != _scatters.end(); ++it) {
      if (**it == ao) {
        _scatters.erase(it);
        break;
      }
    }
  }

  void Analysis::removeAnalysisObject(const Scatter3DPtr& ao) {
    for (vector<shared_ptr<AnalysisObjectPtr> >::iterator it = _scatters.begin();  it != _scatters.end(); ++it) {
      if (**it == ao) {
        _scatters.erase(it);
        break;
      }
    }
  }

  void Analysis::removeAnalysisObject(const MultiweightAOPtr& ao) {
    for (vector<shared_ptr<MultiweightAOPtr> >::iterator it = _analysisobjects.begin();  it != _analysisobjects.end(); ++it) {
      if (**it == ao) {
        _analysisobjects.erase(it);
        break;
      }
    }
  }

}
