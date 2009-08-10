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
#include "Rivet/Analyses/DELPHI_2003_WUD_03_11.hh"
#include "Rivet/Analyses/OPAL_1998_S3780481.hh"
#include "Rivet/Analyses/PDG_Hadron_Multiplicities.hh"
#include "Rivet/Analyses/PDG_Hadron_Multiplicities_Ratios.hh"
//#include "Rivet/Analyses/OPAL_2004_S6132243.hh"
#include "Rivet/Analyses/JADE_OPAL_2000_S4300807.hh"

// HERA
#include "Rivet/Analyses/H1_1994_S2919893.hh"
#include "Rivet/Analyses/H1_1995_S3167097.hh"
#include "Rivet/Analyses/H1_2000_S4129130.hh"
#include "Rivet/Analyses/ZEUS_2001_S4815815.hh"

// Tevatron
#include "Rivet/Analyses/CDF_1988_S1865951.hh"
#include "Rivet/Analyses/CDF_1990_S2089246.hh"
#include "Rivet/Analyses/CDF_1994_S2952106.hh"
#include "Rivet/Analyses/CDF_2000_S4155203.hh"
#include "Rivet/Analyses/CDF_2001_S4751469.hh"
#include "Rivet/Analyses/CDF_2002_S4796047.hh"
#include "Rivet/Analyses/CDF_2004_S5839831.hh"
#include "Rivet/Analyses/CDF_2005_S6080774.hh"
#include "Rivet/Analyses/CDF_2005_S6217184.hh"
#include "Rivet/Analyses/CDF_2006_S6450792.hh"
#include "Rivet/Analyses/CDF_2006_S6653332.hh"
#include "Rivet/Analyses/CDF_2007_S7057202.hh"
#include "Rivet/Analyses/CDF_2008_S7541902.hh"
#include "Rivet/Analyses/CDF_2008_S7540469.hh"
#include "Rivet/Analyses/CDF_2008_NOTE_9351.hh"
#include "Rivet/Analyses/CDF_2008_LEADINGJETS.hh"
#include "Rivet/Analyses/CDF_2008_S7782535.hh"
#include "Rivet/Analyses/CDF_2008_S7828950.hh"
#include "Rivet/Analyses/CDF_2008_S8093652.hh"
#include "Rivet/Analyses/CDF_2008_S8095620.hh"
#include "Rivet/Analyses/CDF_2009_S8233977.hh"
#include "Rivet/Analyses/D0_1996_S3324664.hh"
#include "Rivet/Analyses/D0_2001_S4674421.hh"
#include "Rivet/Analyses/D0_2004_S5992206.hh"
#include "Rivet/Analyses/D0_2006_S6438750.hh"
#include "Rivet/Analyses/D0_2007_S7075677.hh"
#include "Rivet/Analyses/D0_2008_S6879055.hh"
#include "Rivet/Analyses/D0_2008_S7554427.hh"
#include "Rivet/Analyses/D0_2008_S7662670.hh"
#include "Rivet/Analyses/D0_2008_S7719523.hh"
#include "Rivet/Analyses/D0_2008_S7837160.hh"
#include "Rivet/Analyses/D0_2008_S7863608.hh"
#include "Rivet/Analyses/D0_2009_S8202443.hh"
#include "Rivet/Analyses/D0_2009_S8320160.hh"
#include "Rivet/Analyses/D0_2009_S8349509.hh"

// RHIC
#include "Rivet/Analyses/STAR_2006_S6870392.hh"
#include "Rivet/Analyses/STAR_2008_S7993412.hh"
//#include "Rivet/Analyses/STAR_2009_UE_HELEN.hh"

// UA5
#include "Rivet/Analyses/UA5_1986_S1583476.hh" 
#include "Rivet/Analyses/UA5_1988_S1867512.hh"
#include "Rivet/Analyses/UA5_1989_S1926373.hh"
#include "Rivet/Analyses/UA5_1982_S875503.hh"
#include "Rivet/Analyses/UA5_1986_S1583476.hh"

//UA1
#include "Rivet/Analyses/UA1_1990_S2044935.hh"

// MC validation
#include "Rivet/Analyses/MC_TVT1960_ZJETS.hh"
#include "Rivet/Analyses/MC_LHC_LEADINGJETS.hh"
#include "Rivet/Analyses/MC_LHC_ZANALYSIS.hh"
#include "Rivet/Analyses/MC_LHC_DIJET.hh"

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
    fns["DELPHI_2003_WUD_03_11"] = Rivet::DELPHI_2003_WUD_03_11::create;
    fns["OPAL_1998_S3780481"] = Rivet::OPAL_1998_S3780481::create;
    //fns["OPAL_2004_S6132243"] = Rivet::OPAL_2004_S6132243::create;
    fns["JADE_OPAL_2000_S4300807_35GEV"] = Rivet::JADE_OPAL_2000_S4300807_35GEV::create;
    fns["JADE_OPAL_2000_S4300807_44GEV"] = Rivet::JADE_OPAL_2000_S4300807_44GEV::create;
    fns["JADE_OPAL_2000_S4300807_91GEV"] = Rivet::JADE_OPAL_2000_S4300807_91GEV::create;
    fns["JADE_OPAL_2000_S4300807_133GEV"] = Rivet::JADE_OPAL_2000_S4300807_133GEV::create;
    fns["JADE_OPAL_2000_S4300807_161GEV"] = Rivet::JADE_OPAL_2000_S4300807_161GEV::create;
    fns["JADE_OPAL_2000_S4300807_172GEV"] = Rivet::JADE_OPAL_2000_S4300807_172GEV::create;
    fns["JADE_OPAL_2000_S4300807_183GEV"] = Rivet::JADE_OPAL_2000_S4300807_183GEV::create;
    fns["JADE_OPAL_2000_S4300807_189GEV"] = Rivet::JADE_OPAL_2000_S4300807_189GEV::create;

    // HERA
    fns["H1_1994_S2919893"] = Rivet::H1_1994_S2919893::create;
    fns["H1_1995_S3167097"] = Rivet::H1_1995_S3167097::create;
    fns["H1_2000_S4129130"] = Rivet::H1_2000_S4129130::create;
    fns["ZEUS_2001_S4815815"] = Rivet::ZEUS_2001_S4815815::create;

    // Tevatron
    fns["CDF_1988_S1865951"] = Rivet::CDF_1988_S1865951::create;
    fns["CDF_1990_S2089246"] = Rivet::CDF_1990_S2089246::create;
    fns["CDF_1994_S2952106"] = Rivet::CDF_1994_S2952106::create;
    fns["CDF_2000_S4155203"] = Rivet::CDF_2000_S4155203::create;
    fns["CDF_2001_S4751469"] = Rivet::CDF_2001_S4751469::create;
    fns["CDF_2002_S4796047"] = Rivet::CDF_2002_S4796047::create;
    fns["CDF_2004_S5839831"] = Rivet::CDF_2004_S5839831::create;
    fns["CDF_2005_S6080774"] = Rivet::CDF_2005_S6080774::create;
    fns["CDF_2005_S6217184"] = Rivet::CDF_2005_S6217184::create;
    fns["CDF_2006_S6450792"] = Rivet::CDF_2006_S6450792::create;
    fns["CDF_2006_S6653332"] = Rivet::CDF_2006_S6653332::create;
    fns["CDF_2007_S7057202"] = Rivet::CDF_2007_S7057202::create;
    fns["CDF_2008_S7541902"] = Rivet::CDF_2008_S7541902::create;
    fns["CDF_2008_S7540469"] = Rivet::CDF_2008_S7540469::create;
    fns["CDF_2008_S7782535"] = Rivet::CDF_2008_S7782535::create;
    fns["CDF_2008_NOTE_9351"] = Rivet::CDF_2008_NOTE_9351::create;
    fns["CDF_2008_LEADINGJETS"] = Rivet::CDF_2008_LEADINGJETS::create;
    fns["CDF_2008_S7828950"] = Rivet::CDF_2008_S7828950::create;
    fns["CDF_2008_S8093652"] = Rivet::CDF_2008_S8093652::create;
    fns["CDF_2008_S8095620"] = Rivet::CDF_2008_S8095620::create;
    fns["CDF_2009_S8233977"] = Rivet::CDF_2009_S8233977::create;
    fns["D0_2001_S4674421"] = Rivet::D0_2001_S4674421::create;
    fns["D0_1996_S3324664"] = Rivet::D0_1996_S3324664::create;
    fns["D0_2004_S5992206"] = Rivet::D0_2004_S5992206::create;
    fns["D0_2006_S6438750"] = Rivet::D0_2006_S6438750::create;
    fns["D0_2007_S7075677"] = Rivet::D0_2007_S7075677::create;
    fns["D0_2008_S6879055"] = Rivet::D0_2008_S6879055::create;
    fns["D0_2008_S7554427"] = Rivet::D0_2008_S7554427::create;
    fns["D0_2008_S7662670"] = Rivet::D0_2008_S7662670::create;
    fns["D0_2008_S7719523"] = Rivet::D0_2008_S7719523::create;
    fns["D0_2008_S7837160"] = Rivet::D0_2008_S7837160::create;
    fns["D0_2008_S7863608"] = Rivet::D0_2008_S7863608::create;
    fns["D0_2009_S8202443"] = Rivet::D0_2009_S8202443::create;
    fns["D0_2009_S8320160"] = Rivet::D0_2009_S8320160::create;
    fns["D0_2009_S8349509"] = Rivet::D0_2009_S8349509::create;

    // RHIC
    fns["STAR_2006_S6870392"] = Rivet::STAR_2006_S6870392::create;
    fns["STAR_2008_S7993412"] = Rivet::STAR_2008_S7993412::create;
    //fns["STAR_2009_UE_HELEN"] = Rivet::STAR_2009_UE_HELEN::create;

    // UA5
    fns["UA5_1986_S1583476"] = Rivet::UA5_1986_S1583476::create;
    fns["UA5_1988_S1867512"] = Rivet::UA5_1988_S1867512::create;
    fns["UA5_1989_S1926373"] = Rivet::UA5_1989_S1926373::create;
    fns["UA5_1982_S875503"] = Rivet::UA5_1982_S875503::create;
    fns["UA5_1986_S1583476"] = Rivet::UA5_1986_S1583476::create;

   // UA1
    fns["UA1_1990_S2044935"] = Rivet::UA1_1990_S2044935::create;

    // General
    fns["PDG_HADRON_MULTIPLICITIES"] = Rivet::PDG_HADRON_MULTIPLICITIES::create;
    fns["PDG_HADRON_MULTIPLICITIES_RATIOS"] = Rivet::PDG_HADRON_MULTIPLICITIES_RATIOS::create;

    // MC validation
    fns["MC_TVT1960_ZJETS"] = Rivet::MC_TVT1960_ZJETS::create;
    fns["MC_LHC_LEADINGJETS"] = Rivet::MC_LHC_LEADINGJETS::create;
    fns["MC_LHC_ZANALYSIS"] = Rivet::MC_LHC_ZANALYSIS::create;
    fns["MC_LHC_DIJET"] = Rivet::MC_LHC_DIJET::create;

    return fns;
  }

}
