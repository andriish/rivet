// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"

// Example analyses
#include "Rivet/Analyses/ExampleAnalysis.hh"
#include "Rivet/Analyses/ExampleTree.hh"

// LEP
#include "Rivet/Analyses/ALEPH_1991_S2435284.hh"
#include "Rivet/Analyses/ALEPH_1996_S3486095.hh"
#include "Rivet/Analyses/DELPHI_1995_S3137023.hh"
#include "Rivet/Analyses/DELPHI_1996_S3430090.hh"
#include "Rivet/Analyses/DELPHI_2002_069_CONF_603.hh"
#include "Rivet/Analyses/OPAL_1998_S3780481.hh"
#include "Rivet/Analyses/PDG_Hadron_Multiplicities.hh"
#include "Rivet/Analyses/PDG_Hadron_Multiplicities_Ratios.hh"
//#include "Rivet/Analyses/OPAL_2004_S6132243.hh"

// HERA
#include "Rivet/Analyses/H1_1994_S2919893.hh"
#include "Rivet/Analyses/H1_1995_S3167097.hh"
#include "Rivet/Analyses/H1_2000_S4129130.hh"
#include "Rivet/Analyses/ZEUS_2001_S4815815.hh"

// Tevatron
#include "Rivet/Analyses/CDF_1994_S2952106.hh"
#include "Rivet/Analyses/CDF_2001_S4751469.hh"
#include "Rivet/Analyses/CDF_2004_S5839831.hh"
#include "Rivet/Analyses/CDF_2005_S6217184.hh"
#include "Rivet/Analyses/CDF_2006_S6653332.hh"
#include "Rivet/Analyses/CDF_2007_S7057202.hh"
#include "Rivet/Analyses/CDF_2008_S7541902.hh"
#include "Rivet/Analyses/CDF_2008_MINBIAS.hh"
#include "Rivet/Analyses/CDF_2008_DRELLYAN.hh"
#include "Rivet/Analyses/D0_2001_S4674421.hh"
#include "Rivet/Analyses/D0_2004_S5992206.hh"
#include "Rivet/Analyses/D0_2008_S6879055.hh"


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
    fns["ALEPH_1996_S3486095"] = Rivet::ALEPH_1996_S3486095::create;
    fns["DELPHI_1995_S3137023"] = Rivet::DELPHI_1995_S3137023::create;
    fns["DELPHI_1996_S3430090"] = Rivet::DELPHI_1996_S3430090::create;
    fns["DELPHI_2002_069_CONF_603"] = Rivet::DELPHI_2002_069_CONF_603::create;
    fns["OPAL_1998_S3780481"] = Rivet::OPAL_1998_S3780481::create;
    //fns["OPAL_2004_S6132243"] = Rivet::OPAL_2004_S6132243::create;

    // HERA
    fns["H1_1994_S2919893"] = Rivet::H1_1994_S2919893::create;
    fns["H1_1995_S3167097"] = Rivet::H1_1995_S3167097::create;
    fns["H1_2000_S4129130"] = Rivet::H1_2000_S4129130::create;
    fns["ZEUS_2001_S4815815"] = Rivet::ZEUS_2001_S4815815::create;

    // Tevatron
    fns["CDF_1994_S2952106"] = Rivet::CDF_1994_S2952106::create;
    fns["CDF_2001_S4751469"] = Rivet::CDF_2001_S4751469::create;
    fns["CDF_2004_S5839831"] = Rivet::CDF_2004_S5839831::create;
    fns["CDF_2005_S6217184"] = Rivet::CDF_2005_S6217184::create;
    fns["CDF_2006_S6653332"] = Rivet::CDF_2006_S6653332::create;
    fns["CDF_2007_S7057202"] = Rivet::CDF_2007_S7057202::create;
    fns["CDF_2008_S7541902"] = Rivet::CDF_2008_S7541902::create;
    fns["CDF_2008_MINBIAS"] = Rivet::CDF_2008_MINBIAS::create;
    fns["CDF_2008_DRELLYAN"] = Rivet::CDF_2008_DRELLYAN::create;
    fns["D0_2001_S4674421"] = Rivet::D0_2001_S4674421::create;
    fns["D0_2004_S5992206"] = Rivet::D0_2004_S5992206::create;
    fns["D0_2008_S6879055"] = Rivet::D0_2008_S6879055::create;

    // General
    fns["PDG_HADRON_MULTIPLICITIES"] = Rivet::PDG_Hadron_Multiplicities::create;
    fns["PDG_HADRON_MULTIPLICITIES_RATIOS"] = Rivet::PDG_Hadron_Multiplicities_Ratios::create;

    return fns;
  }

}
