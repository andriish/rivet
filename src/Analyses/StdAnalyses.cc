// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"

// Concrete analyses
#include "Rivet/Analyses/TestAnalysis.hh"
////#include "Rivet/Analyses/EurPhys40C287.hh"
#include "Rivet/Analyses/ExampleTree.hh"
//#include "Rivet/Analyses/HepEx0112029.hh"
#include "Rivet/Analyses/HepEx0409040.hh"
#include "Rivet/Analyses/HepEx0505013.hh"
//#include "Rivet/Analyses/HepEx9506012.hh"
#include "Rivet/Analyses/HepEx0701051.hh"
#include "Rivet/Analyses/PL273B181.hh"
#include "Rivet/Analyses/PRD65092002.hh"
#include "Rivet/Analyses/ZPhys73C11.hh"
#include "Rivet/Analyses/HepEx0605099.hh"
//#include "Rivet/Analyses/MyAnalysis.hh"


extern "C" {

  /// Get a map of factory functions for making analyses, keyed 
  /// by upper-case versions of the analyses' names. (Upper-case 
  /// keys are used since they can be uniquely obtained from any 
  /// mixed-case key, which allows a bit of UI flexibility.)
  AnalysisBuilders getAnalysisBuilders() {
    AnalysisBuilders fns;
    fns["TEST"] = Rivet::TestAnalysis::create;
    ////fns["EURPHYS40C287"] = Rivet::EurPhys40C287::create;
    fns["EXAMPLETREE"] = Rivet::ExampleTree::create;
    //fns["HEPEX0112029"] = Rivet::HepEx0112029::create;
    fns["HEPEX0409040"] = Rivet::HepEx0409040::create;
    fns["HEPEX0505013"] = Rivet::HepEx0505013::create;
    //fns["HEPEX9506012"] = Rivet::HepEx9506012::create;
    fns["HEPEX0701051"] = Rivet::HepEx0701051::create;
    fns["PL273B181"] = Rivet::PL273B181::create;
    fns["PRD65092002"] = Rivet::PRD65092002::create;
    fns["ZPHYS73C11"] = Rivet::ZPhys73C11::create;
    fns["HEPEX0605099"] = Rivet::HepEx0605099::create;
    //fns["MYANALYSIS"] = Rivet::MyAnalysis::create;
    return fns;
  }

}
