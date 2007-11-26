#ifndef RIVETGUN_AVAILABLEGENERATORS_HH 
#define RIVETGUN_AVAILABLEGENERATORS_HH 

#include "AGILe/Generator.hh"

#ifdef HAVE_GENERATOR_FPYTHIA
#include "AGILe/FPythia/FPythia.hh"
#endif

#ifdef HAVE_GENERATOR_FHERWIG
#include "AGILe/FHerwig/FHerwig.hh"
#endif

#ifdef HAVE_GENERATOR_FALPHER
#include "AGILe/AlpGen/AlpGenFHerwig.hh"
#endif

#ifdef HAVE_GENERATOR_FALPPYT
#include "AGILe/AlpGen/AlpGenFPythia.hh"
#endif

#ifdef HAVE_GENERATOR_CHARYBDISHER
#include "AGILe/Charybdis/CharybdisFHerwig.hh"
#endif

#ifdef HAVE_GENERATOR_CHARYBDISPYT
#include "AGILe/Charybdis/CharybdisFPythia.hh"
#endif

#ifdef HAVE_GENERATOR_SHERPA
#include "AGILe/Sherpa/Sherpa.hh"
#endif

#ifdef HAVE_GENERATOR_CCHERWIG
#include "AGILe/CCHerwig/CCHerwig.hh"
#endif

#ifdef HAVE_GENERATOR_CCPYTHIA
#include "AGILe/CCPythia/CCPythia.hh"
#endif


namespace Rivet {

  inline vector<string> getAvailableGenNames() {
    vector<string> genNames;
    #ifdef HAVE_GENERATOR_FPYTHIA
    genNames.push_back("FPythia");
    #endif
    #ifdef HAVE_GENERATOR_FHERWIG
    genNames.push_back("FHerwig");
    #endif
    #ifdef HAVE_GENERATOR_FALPHER
    genNames.push_back("AlpGenFHerwig");
    #endif
    #ifdef HAVE_GENERATOR_FALPPYT
    genNames.push_back("AlpGenFPythia");
    #endif
    #ifdef HAVE_GENERATOR_CHARYBDISHER
    genNames.push_back("CharybdisFHerwig");
    #endif
    #ifdef HAVE_GENERATOR_CHARYBDISPYT
    genNames.push_back("CharybdisFPythia");
    #endif
    #ifdef HAVE_GENERATOR_SHERPA
    genNames.push_back("Sherpa");
    #endif
    #ifdef HAVE_GENERATOR_CCHERWIG
    genNames.push_back("CCHerwig");
    #endif
    #ifdef HAVE_GENERATOR_CCPYTHIA
    genNames.push_back("CCPythia");
    #endif

    return genNames;
  }


  inline AGILe::Generator* makeNewGenerator(string name) {
    #ifdef HAVE_GENERATOR_FHERWIG
    if (name == "FHerwig") return new AGILe::FHerwig();
    #ifdef HAVE_GENERATOR_FALPHER
    if (name == "AlpGenFHerwig") return new AGILe::AlpGenFHerwig();
    #endif
    #ifdef HAVE_GENERATOR_CHARYBDISHER
    if (name == "CharybdisFHerwig") return new AGILe::CharybdisFHerwig();
    #endif
    #endif

    #ifdef HAVE_GENERATOR_FPYTHIA
    if (name == "FPythia") return new AGILe::FPythia();
    #ifdef HAVE_GENERATOR_FALPPYT
    if (name == "AlpGenFPythia") return new AGILe::AlpGenFPythia();
    #endif
    #ifdef HAVE_GENERATOR_CHARYBDISPYT
    if (name == "CharybdisFPythia") return new AGILe::CharybdisFPythia();
    #endif
    #endif

    #ifdef HAVE_GENERATOR_SHERPA
    if (name == "Sherpa") return new AGILe::Sherpa();
    #endif

    #ifdef HAVE_GENERATOR_CCHERWIG
    if (name == "CCHerwig") return new AGILe::CCHerwig();
    #endif

    #ifdef HAVE_GENERATOR_CCPYTHIA
    if (name == "CCPythia") return new AGILe::CCPythia();
    #endif

    return 0;
  }


}

#endif
