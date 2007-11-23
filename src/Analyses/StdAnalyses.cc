// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"

// Concrete analyses
#include "Rivet/Analyses/TestAnalysis.hh"
////#include "Rivet/Analyses/S6132243.hh"
#include "Rivet/Analyses/ExampleTree.hh"
//#include "Rivet/Analyses/S4815815.hh"
#include "Rivet/Analyses/S5992206.hh"
#include "Rivet/Analyses/S6217184.hh"
//#include "Rivet/Analyses/S3167097.hh"
#include "Rivet/Analyses/S7057202.hh"
#include "Rivet/Analyses/S2435284.hh"
#include "Rivet/Analyses/S4751469.hh"
#include "Rivet/Analyses/S3430090.hh"
#include "Rivet/Analyses/S6653332.hh"
//#include "Rivet/Analyses/MyAnalysis.hh"


extern "C" {

  /// Get a map of factory functions for making analyses, keyed 
  /// by upper-case versions of the analyses' names. (Upper-case 
  /// keys are used since they can be uniquely obtained from any 
  /// mixed-case key, which allows a bit of UI flexibility.)
  AnalysisBuilders getAnalysisBuilders() {
    AnalysisBuilders fns;
    fns["TEST"] = Rivet::TestAnalysis::create;
    ////fns["EURPHYS40C287"] = Rivet::S6132243::create;
    fns["EXAMPLETREE"] = Rivet::ExampleTree::create;
    //fns["HEPEX0112029"] = Rivet::S4815815::create;
    fns["HEPEX0409040"] = Rivet::S5992206::create;
    fns["HEPEX0505013"] = Rivet::S6217184::create;
    //fns["HEPEX9506012"] = Rivet::S3167097::create;
    fns["HEPEX0701051"] = Rivet::S7057202::create;
    fns["S2435284"] = Rivet::S2435284::create;
    fns["S4751469"] = Rivet::S4751469::create;
    fns["ZPHYS73C11"] = Rivet::S3430090::create;
    fns["HEPEX0605099"] = Rivet::S6653332::create;
    //fns["MYANALYSIS"] = Rivet::MyAnalysis::create;
    return fns;
  }

}
