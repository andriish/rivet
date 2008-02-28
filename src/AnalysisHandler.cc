// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/RivetAIDA.hh"
#include "AIDA/IManagedObject.h"

using namespace AIDA;


namespace Rivet {


  AnalysisHandler::AnalysisHandler(string basefilename, HistoFormat storetype)
    : _nRun(0), _iRun(0), _numEvents(0), _sumOfWeights(0.0) {
    _theAnalysisFactory = AIDA_createAnalysisFactory();
    setupFactories(basefilename, storetype);
  }


  AnalysisHandler::AnalysisHandler(IAnalysisFactory& afac, string basefilename, HistoFormat storetype)
    : _nRun(0), _iRun(0), _numEvents(0), _sumOfWeights(0.0), 
      _theAnalysisFactory(&afac) {
    setupFactories(basefilename, storetype);
  }


  Log& AnalysisHandler::getLog() { 
    return Log::getLog("Rivet.Analysis.Handler");
  }


  void AnalysisHandler::init(int i, int N) {
    _nRun = N;
    _iRun = i;
    _numEvents = 0;
    _sumOfWeights = 0.0;
    for (set<Analysis*>::iterator a = _analyses.begin(); a != _analyses.end(); ++a) {
      getLog() << Log::DEBUG << "Initialising analysis: " << (*a)->getName() << endl;
      (*a)->init();
      getLog() << Log::DEBUG << "Checking consistency of analysis: " << (*a)->getName() << endl;
      (*a)->checkConsistency();
      getLog() << Log::DEBUG << "Done initialising analysis: " << (*a)->getName() << endl;
    }
  }
  

  void AnalysisHandler::analyze(const GenEvent& geneve) {
    Event event(geneve);
    _numEvents++;
    _sumOfWeights += event.weight();
    for (set<Analysis*>::iterator a = _analyses.begin(); a != _analyses.end(); ++a) {
      getLog() << Log::DEBUG << "About to run analysis " << (*a)->getName() << endl;
      (*a)->analyze(event);
      getLog() << Log::DEBUG << "Finished running analysis " << (*a)->getName() << endl;
    }
  }


  void AnalysisHandler::finalize() {
    Log& log = getLog();
    log << Log::INFO << "Finalising analysis" << endl;
    for (set<Analysis*>::iterator a = _analyses.begin(); a != _analyses.end(); ++a) {
      (*a)->finalize();
    }

    // Change AIDA histos into data point sets
    log << Log::INFO << "Normalising the AIDA tree" << endl;
    assert(_theTree != 0);
    normalizeTree(tree());
    tree().commit();

    // Delete analyses
    for (set<Analysis*>::iterator a = _analyses.begin(); a != _analyses.end(); ++a) {
      delete *a;
    }
    _analyses.clear();
    AnalysisLoader::closeAnalysisBuilders();
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname) {
    Analysis* analysis = AnalysisLoader::getAnalysis(analysisname);
    if (analysis) { // < Check for null analysis.
      analysis->_theHandler = this;
      _analyses.insert(analysis);
    }
    return *this;
  }


  void AnalysisHandler::setupFactories(string basefilename, HistoFormat storetype) {
    string filename(basefilename), storetypestr("");
    if (storetype == AIDAML) {
      filename += ".aida";
      storetypestr = "xml";
    } else if (storetype == FLAT) {
      filename += ".data";
      storetypestr = "flat";
    } else if (storetype == ROOT) {
      filename += ".root";
      storetypestr = "root";
    }
    _theTree = _theAnalysisFactory->createTreeFactory()->create(filename, storetypestr, false, true);
    _theHistogramFactory = _theAnalysisFactory->createHistogramFactory(tree());
    _theDataPointSetFactory = _theAnalysisFactory->createDataPointSetFactory(tree());
  }


  void AnalysisHandler::normalizeTree(ITree& tree) {
    Log& log = getLog();
    const vector<string> paths = tree.listObjectNames("/", true); // args set recursive listing
    log << Log::TRACE << "Number of objects in AIDA tree = " << paths.size() << endl;
    const string tmpdir = "/RivetNormalizeTmp";
    tree.mkdir(tmpdir);
    for (vector<string>::const_iterator path = paths.begin(); path != paths.end(); ++path) {
      
      IManagedObject* hobj = tree.find(*path);
      if (hobj) {
        IHistogram1D* histo = dynamic_cast<IHistogram1D*>(hobj);
        IProfile1D* prof = dynamic_cast<IProfile1D*>(hobj);
        // If it's a normal histo:
        if (histo) {
          log << Log::TRACE << "Converting histo " << *path << " to DPS" << endl;
          tree.mv(*path, tmpdir);
          const size_t lastslash = path->find_last_of("/");
          const string basename = path->substr(lastslash+1, path->length() - (lastslash+1));
          const string tmppath = tmpdir + "/" + basename;
          IHistogram1D* tmphisto = dynamic_cast<IHistogram1D*>(tree.find(tmppath));
          if (tmphisto) {
            //getLog() << Log::TRACE << "Temp histo " << tmppath << " exists" << endl;
            datapointsetFactory().create(*path, *tmphisto);
          }
          tree.rm(tmppath);
        }
        // If it's a profile histo:
        else if (prof) {
          log << Log::TRACE << "Converting profile histo " << *path << " to DPS" << endl;
          tree.mv(*path, tmpdir);
          const size_t lastslash = path->find_last_of("/");
          const string basename = path->substr(lastslash+1, path->length() - (lastslash+1));
          const string tmppath = tmpdir + "/" + basename;
          IProfile1D* tmpprof = dynamic_cast<IProfile1D*>(tree.find(tmppath));
          if (tmpprof) {
            //getLog() << Log::TRACE << "Temp profile histo " << tmppath << " exists" << endl;
            datapointsetFactory().create(*path, *tmpprof);
          }
          tree.rm(tmppath);
        }

      }
      
      
    }
    tree.rmdir(tmpdir);
  }


}
