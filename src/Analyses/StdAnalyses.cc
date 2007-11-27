// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"

// Example analyses
#include "Rivet/Analyses/ExampleAnalysis.hh"
#include "Rivet/Analyses/ExampleTree.hh"

// LEP
#include "Rivet/Analyses/ALEPH_1991_S2435284.hh"
#include "Rivet/Analyses/DELPHI_1996_S3430090.hh"
#include "Rivet/Analyses/OPAL_2004_S6132243.hh"

// HERA
#include "Rivet/Analyses/H1_1995_S3167097.hh"
#include "Rivet/Analyses/ZEUS_2001_S4815815.hh"

// Tevatron
#include "Rivet/Analyses/CDF_2001_S4751469.hh"
#include "Rivet/Analyses/CDF_2005_S6217184.hh"
#include "Rivet/Analyses/CDF_2006_S6653332.hh"
#include "Rivet/Analyses/CDF_2007_S7057202.hh"
//#include "Rivet/Analyses/D0_2001_S4674421.hh"
#include "Rivet/Analyses/D0_2004_S5992206.hh"



extern "C" {

  /// Get a map of factory functions for making analyses, keyed 
  /// by upper-case versions of the analyses' names. (Upper-case 
  /// keys are used since they can be uniquely obtained from any 
  /// mixed-case key, which allows a bit of UI flexibility.)
  AnalysisBuilders getAnalysisBuilders() {
    AnalysisBuilders fns;
    fns["EXAMPLE"] = Rivet::ExampleAnalysis::create;
    fns["EXAMPLETREE"] = Rivet::ExampleTree::create;

    // LEP
    fns["ALEPH_1991_S2435284"] = Rivet::ALEPH_1991_S2435284::create;
    fns["DELPHI_1996_S3430090"] = Rivet::DELPHI_1996_S3430090::create;
    fns["OPAL_2004_S6132243"] = Rivet::OPAL_2004_S6132243::create;

    // HERA
    fns["H1_1995_S3167097"] = Rivet::H1_1995_S3167097::create;
    fns["ZEUS_2001_S4815815"] = Rivet::ZEUS_2001_S4815815::create;

    // Tevatron
    fns["CDF_2001_S4751469"] = Rivet::CDF_2001_S4751469::create;
    fns["CDF_2005_S6217184"] = Rivet::CDF_2005_S6217184::create;
    fns["CDF_2006_S6653332"] = Rivet::CDF_2006_S6653332::create;
    fns["CDF_2007_S7057202"] = Rivet::CDF_2007_S7057202::create;
    //fns["D0_2001_S4674421"] = Rivet::D0_2001_S4674421::create;
    fns["D0_2004_S5992206"] = Rivet::D0_2004_S5992206::create;

    return fns;
  }

}
