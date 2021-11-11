// main89.cc is a part of the PYTHIA event generator.
// Copyright (C) 2021 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Stefan Prestel <stefan.prestel@thep.lu.se>.

// Keywords: matching; merging; leading order; NLO; powheg; madgraph; aMC@NLO;
// CKKW-L; UMEPS; NL3; UNLOPS; FxFx; MLM; userhooks; LHE file; hepmc;

// This program illustrates how to do run PYTHIA with LHEF input, allowing a
// sample-by-sample generation of
// a) Non-matched/non-merged events
// b) MLM jet-matched events (kT-MLM, shower-kT, FxFx)
// c) CKKW-L and UMEPS-merged events
// d) UNLOPS NLO merged events
// see the respective sections in the online manual for details.
//
// An example command is
//     ./main89 main89ckkwl.cmnd hepmcout89.dat
// where main89.cmnd supplies the commands and hepmcout89.dat is the
// output file. This example requires HepMC2 or HepMC3.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"


#include "Pythia8Plugins/Pythia8Rivet.h"
#include "main93.h"


// Include UserHooks for Jet Matching.
#include "Pythia8Plugins/CombineMatchingInput.h"
// Include UserHooks for randomly choosing between integrated and
// non-integrated treatment for unitarised merging.
#include "Pythia8Plugins/aMCatNLOHooks.h"

using namespace Pythia8;

//==========================================================================

// Example main programm to illustrate merging.

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 2) {
    cerr << " Unexpected number of command-line arguments ("<<argc<<"). \n"
         << " You are expected to provide argument" << endl
         << " 1. Input file for settings" << endl
         << " Program stopped. " << endl;
    return 1;
  }
  
  string out = "";


  Pythia pythia;
  
  // UserHooks wrapper
  auto userHooksWrapper = make_shared<UserHooksWrapper>();
  userHooksWrapper->additionalSettings(&pythia.settings);
  pythia.setUserHooksPtr(userHooksWrapper);
  pythia.settings.addFlag("Main:writeHepMC",false);
  pythia.settings.addFlag("Main:writeRoot",false);
  pythia.settings.addFlag("Main:runRivet",false);
  pythia.settings.addFlag("Main:rivetIgnoreBeams",false);
  pythia.settings.addFlag("Main:outputLog",false);
  pythia.settings.addWVec("Main:analyses",vector<string>());
  pythia.settings.addWVec("Main:preload",vector<string>());
  pythia.settings.addWord("Main:rivetRunName","");

  // New setting to allow processing of multiple input LHEFs.
  pythia.settings.addMode("LHEFInputs:nSubruns",0,true,false,0,100);

  // Input parameters:
  pythia.readFile(argv[1],0);

  // Deactivate AUX_ weight output
  pythia.readString("Weights:suppressAUX = on");
  
  const bool hepmc = pythia.flag("Main:writeHepMC");
  const bool root = pythia.flag("Main:writeRoot");
  const bool runRivet = pythia.flag("Main:runRivet");
  const bool ignoreBeams = pythia.flag("Main:rivetIgnoreBeams");
  const bool doLog = pythia.flag("Main:outputLog");
  const string rivetrName = pythia.settings.word("Main:rivetRunName");
  const vector<string> rAnalyses = pythia.settings.wvec("Main:analyses");
  const vector<string> rPreload = pythia.settings.wvec("Main:preload");
  int nError = pythia.mode("Main:timesAllowErrors");
  
  // Rivet initialization.
  Pythia8Rivet rivet(pythia,(out == "" ? "Rivet.yoda" : out + ".yoda"));
  for(int i = 0, N = rAnalyses.size(); i < N; ++i){
    string analysis = rAnalyses[i];
    size_t pos = analysis.find(":");
    // Simple case, no analysis parameters.
    if(pos == string::npos)
      rivet.addAnalysis(analysis);
    else {
      string an = analysis.substr(0,pos);
      analysis.erase(0, pos + 1);
      vector<string> pKeys;
      vector<string> pVals;
      while (analysis.find("->") != string::npos) {
        pos = analysis.find(":");
        string par = analysis.substr(0,pos);
        size_t pos2 = par.find("->");
        if (pos2 == string::npos){
           cout << "Error: malformed parameter " << par << endl;
        }
        string pKey = par.substr(0,pos2);
        string pVal = par.substr(pos2+2,par.length());
        pKeys.push_back(pKey);
        pVals.push_back(pVal);
        analysis.erase(0,par.length()+1);
      }
      for (int j = 0, N = pKeys.size(); j < N; ++j)
        an += ":"+pKeys[j]+"="+pVals[j];
      rivet.addAnalysis(an);
    }
  }
  for(int i = 0, N = rPreload.size(); i < N; ++i)
    rivet.addPreload(rPreload[i]);
  rivet.addRunName(rivetrName);
  
  // Check if jet matching should be applied.
  bool doMatch   = pythia.settings.flag("JetMatching:merge");

  // Check if internal merging should be applied.
  bool doMerge   = !(pythia.settings.word("Merging:Process").compare("void")
    == 0);

  // Currently, only one scheme at a time is allowed.
  if (doMatch && doMerge) {
    cerr << " Jet matching and merging cannot be used simultaneously.\n"
         << " Program stopped.";
  }

  // Get number of subruns.
  int nMerge = pythia.mode("LHEFInputs:nSubruns");
  bool doMatchMerge = true;
  if (nMerge == 0) { nMerge = 1; doMatchMerge = false; }

  // Number of events. Negative numbers mean all events in the LHEF will be
  // used.
  long nEvent = pythia.settings.mode("Main:numberOfEvents");
  if (nEvent < 1) nEvent = 1000;

  // For jet matching, initialise the respective user hooks code.
  //shared_ptr<UserHooks> matching;

  // Allow to set the number of addtional partons dynamically.
  shared_ptr<amcnlo_unitarised_interface> setting;
  if ( doMerge ) {
    // Store merging scheme.
    int scheme = ( pythia.settings.flag("Merging:doUMEPSTree")
                || pythia.settings.flag("Merging:doUMEPSSubt")) ?
                1 :
                 ( ( pythia.settings.flag("Merging:doUNLOPSTree")
                || pythia.settings.flag("Merging:doUNLOPSSubt")
                || pythia.settings.flag("Merging:doUNLOPSLoop")
                || pythia.settings.flag("Merging:doUNLOPSSubtNLO")) ?
                2 :
                0 );
    setting = make_shared<amcnlo_unitarised_interface>(scheme);
    pythia.setUserHooksPtr(setting);
  }

  // For jet matching, initialise the respective user hooks code.
  CombineMatchingInput combined;
  if (doMatch) combined.setHook(pythia);

  vector<double> xss;

  // Allow usage also for non-matched configuration.
  if(!doMatchMerge) {
    // Loop over subruns with varying number of jets.
    for (int iMerge = 0; iMerge < nMerge; ++iMerge) {
      // Read in file for current subrun and initialize.
      pythia.readFile(argv[1], iMerge);
      // Initialise.
      pythia.init();
      // Start generation loop
      while( pythia.info.nSelected() < nEvent ){
        // Generate next event
        if( !pythia.next() ) {
          if ( pythia.info.atEndOfFile() ) break;
          else continue;
        }
      } // end loop over events to generate.
      // print cross section, errors
      pythia.stat();
      xss.push_back(pythia.info.sigmaGen());
    }
    pythia.info.weightContainerPtr->clearTotal();
  }

  // Allow abort of run if many errors.
  int  nAbort  = pythia.mode("Main:timesAllowErrors");
  int  iAbort  = 0;
  bool doAbort = false;

  cout << endl << endl << endl;
  cout << "Start generating events" << endl;

  // Loop over subruns with varying number of jets.
  for (int iMerge = 0; iMerge < nMerge; ++iMerge) {

    // Read in name of LHE file for current subrun and initialize.
    pythia.readFile(argv[1], iMerge);

    // Initialise.
    pythia.init();

    // Get the inclusive x-section by summing over all process x-sections.
    double xs = 0.;
    for (int i=0; i < pythia.info.nProcessesLHEF(); ++i)
      xs += pythia.info.sigmaLHEF(i);

    if (!doMatchMerge) xs = xss[iMerge];

    // Start generation loop
    while( pythia.info.nSelected() < nEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if ( pythia.info.atEndOfFile() ) break;
        else if (++iAbort > nAbort) {doAbort = true; break;}
        else continue;
      }

      // Get event weight(s).
      // Additional weight due to random choice of reclustered/non-reclustered
      // treatment. Also contains additional sign for subtractive samples.
      double evtweight = pythia.info.weightValueByIndex();

      // Do not print zero-weight events.
      if ( evtweight == 0. ) continue;

      // Do not print broken / empty events
      if (pythia.event.size() < 3) continue;

      // Work with weighted (LHA strategy=-4) events.
      double norm = 1.;
      if (abs(pythia.info.lhaStrategy()) == 4)
        norm = 1. / double(1e9*nEvent);
      // Work with unweighted events.
      else
        norm = xs / double(1e9*nEvent);

      pythia.info.weightContainerPtr->accumulateXsec(norm);

    if (runRivet) {
      rivet();
    }


    } // end loop over events to generate.

    if (doAbort) break;

    // Print cross section, errors.
    pythia.stat();

    // Get cross section statistics for sample.
    double sigmaSample = pythia.info.weightContainerPtr->getSampleXsec()[0];
    double errorSample = pythia.info.weightContainerPtr->getSampleXsecErr()[0];

    cout << endl << " Contribution of sample " << iMerge
         << " to the inclusive cross section : "
         << scientific << setprecision(8)
         << sigmaSample << "  +-  " << errorSample  << endl;
  }
  cout << endl << endl << endl;

  // Get cross section statistics for total run.
  double sigmaTotal = pythia.info.weightContainerPtr->getTotalXsec()[0];
  double errorTotal = pythia.info.weightContainerPtr->getSampleXsecErr()[0];
  if (doAbort)
    cout << " Run was not completed owing to too many aborted events" << endl;
  else
    cout << "Inclusive cross section: " << scientific << setprecision(8)
         << sigmaTotal << "  +-  " << errorTotal << " mb " << endl;
  cout << endl << endl << endl;

  // Done
  return 0;

}
