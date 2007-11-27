// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"

// Concrete analyses
#include "Rivet/Analyses/ExampleAnalysis.hh"
#include "Rivet/Analyses/ExampleTree.hh"

#include "Rivet/Analyses/S2435284.hh"
//#include "Rivet/Analyses/S3167097.hh"
#include "Rivet/Analyses/S3430090.hh"
#include "Rivet/Analyses/S4751469.hh"
//#include "Rivet/Analyses/S4815815.hh"
#include "Rivet/Analyses/S5992206.hh"
//#include "Rivet/Analyses/S6132243.hh"
#include "Rivet/Analyses/S6217184.hh"
#include "Rivet/Analyses/S6653332.hh"
#include "Rivet/Analyses/S7057202.hh"


extern "C" {

  /// Get a map of factory functions for making analyses, keyed 
  /// by upper-case versions of the analyses' names. (Upper-case 
  /// keys are used since they can be uniquely obtained from any 
  /// mixed-case key, which allows a bit of UI flexibility.)
  AnalysisBuilders getAnalysisBuilders() {
    AnalysisBuilders fns;
    fns["EXAMPLE"] = Rivet::ExampleAnalysis::create;
    fns["EXAMPLETREE"] = Rivet::ExampleTree::create;
    fns["S2435284"] = Rivet::S2435284::create;
    //fns["S3167097"] = Rivet::S3167097::create;
    fns["S3430090"] = Rivet::S3430090::create;
    fns["S4751469"] = Rivet::S4751469::create;
    //fns["S4815815"] = Rivet::S4815815::create;
    fns["S5992206"] = Rivet::S5992206::create;
    //fns["S6132243"] = Rivet::S6132243::create;
    fns["S6217184"] = Rivet::S6217184::create;
    fns["S6653332"] = Rivet::S6653332::create;
    fns["S7057202"] = Rivet::S7057202::create;
    return fns;
  }

}
