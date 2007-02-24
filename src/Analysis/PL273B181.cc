// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/PL273B181.hh"
#include "HepPDT/ParticleID.hh"
#include "AIDA/IAxis.h"

// includes for chi squared calculations:
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <cctype>


using namespace Rivet;
using namespace HepMC;
using namespace std;


PL273B181::~PL273B181() {}

//------------------------------------------------
// Functions for chi squared calculations:

// Structure for storing the data
struct histogramData {
  double bin;
  double yExp;
  double yExpErr;  
  };
 
// Function to split the lines of the data files
vector<string> split(const string& s)
{
  vector<string> ret;
  typedef string::size_type string_size;
  string_size i = 0;

  while (i != s.size()) {
    while (i != s.size() && isspace(s[i]))
      ++i;
      string_size j = i;
      while (j != s.size() && !isspace(s[j]))
	++j;
	if (i != j) {
	  ret.push_back(s.substr(i, j - i));
	  i = j;
	}
  }
  return ret;
}

// Calculate Chi Squared for a given set of values/error
double chi (double yexp, double yerr, double ysim) {
  if (yerr == 0) {
    return ((ysim - yexp)*(ysim-yexp)); 
  }
  else {
    return  ((ysim - yexp)*(ysim-yexp)/(yerr*yerr));
  }
}

// Template function to convert to double
template <class T>
double conv (T value) {
  double doubleOut;   stringstream valueIn;
  valueIn.clear();
  valueIn << value;
  valueIn >> doubleOut;
  return doubleOut;
  }

// Main function that will calculate chi squared 
// from two given histograms (in files)
double ChiError(const char* simData, const char* expData)
{
  string s, test;
  vector<string> store;
  vector<double> ySim;
  histogramData hD;
  vector<histogramData> mExp;
  double chitot = 0;
  ySim.clear();
  ifstream expFile(expData);
  ifstream simFile(simData);

  //read in experimental data for comparison with all other files
  while (getline(expFile, s)) {                           
    store.clear();
    store = split(s);
    test = store[0];
    if (isalnum(test[0])) {
      hD.bin     = conv(store[0]);
      hD.yExp    = conv(store[1]);
      hD.yExpErr = (conv(store[2])+conv(store[3]));
      mExp.push_back(hD);
    }
  } 
     
  //read in the simulated data file and store it for comparison       
  while (getline(simFile, s)) {
    store.clear();
    store = split(s);
    test = store[0];
    if (isalnum(test[0])) {
      ySim.push_back(conv(store[1]));
    }
  }

  //add up the chi squared values
  for (vector<double>::size_type i = 0; i != ySim.size(); ++i) {
    chitot += chi(mExp[i].yExp, mExp[i].yExpErr, ySim[i]);          
 
 cout << "  yExp:   = " << mExp[i].yExp
      << "  yErr:    = " << mExp[i].yExpErr     
      << "  ySim:   = " << ySim[i] 
      << endl;

  }
return chitot;
}

// ------------------------------------------------


// Book histograms
void PL273B181::init() {

  // Book histograms
  AIDA::IHistogramFactory* hf = histogramFactory();
  tree()->mkdir("/Test/");
  histChTot_       = hf->createHistogram1D("/TotalChMult","Total charged multiplicity", 25, 1.0, 51.0);
  histSphericity_  = hf->createHistogram1D("/Sphericity", "Event Shape: Sphericity", 8, 0.0, 0.70);
  histAplanarity_  = hf->createHistogram1D("/Aplanarity", "Event Shape: APlanarity", 10, 0.0, 0.09);
  histPlanarity_   = hf->createHistogram1D("/Planarity",  "Event Shape: Planarity", 16, 0.0, 0.70);
}

// Do the analysis
void PL273B181::analyze(const Event & event) {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);
  log << LogPriority::DEBUG << "Starting analyzing" << endlog;

  // Analyse and print some info
  const Multiplicity& m = event.addProjection(mult);
  //log << LogPriority::INFO << "Total charged multiplicity    = " 
  //    << m.totalChargedMultiplicity() << endlog;

  //Analyse the event shape info
  const Sphericity& s = event.addProjection(spher);
  log << LogPriority::INFO << "Sphericity    = " 
      << s.eventSphericity() << endlog;
  log << LogPriority::INFO << "Aplanarity    = " 
      << s.eventAplanarity() << endlog;
  log << LogPriority::INFO << "Planarity     = " 
      << s.eventPlanarity() << endlog;

  // Fill histograms here, and scale them later
  histChTot_->fill(m.totalChargedMultiplicity(), 1.0);
  histSphericity_->fill(s.eventSphericity(), 1.0);
  histPlanarity_->fill(s.eventPlanarity(), 1.0);
  histAplanarity_->fill(s.eventAplanarity(), 1.0);
  
  // Finished...
  log << LogPriority::DEBUG << "Finished analyzing" << endlog;
}

// Finalize
void PL273B181::finalize() { 
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);

  double area = 0;
  for (int i=0; i < histSphericity_->axis().bins(); ++i) {
    area += histSphericity_->binHeight(i) * histSphericity_->axis().binWidth(i); 
  }
  log << LogPriority::INFO << "Area under histogram: " << area << endlog;


  // Normalize the histogram areas to 1
  histChTot_->scale(1/histChTot_->sumBinHeights() );
  //histSphericity_->scale(1/histSphericity_->sumBinHeights() );
  histSphericity_->scale(1/area);
  histPlanarity_->scale(1/histPlanarity_->sumBinHeights() ); 
  histAplanarity_->scale(1/histAplanarity_->sumBinHeights() );

  log << LogPriority::INFO << "Sum of sph bin heights after normalization: " 
      << histSphericity_->sumBinHeights() << endlog;
}


RivetInfo PL273B181::getInfo() const {
  return Analysis::getInfo() + fsproj.getInfo() +
    mult.getInfo() + spher.getInfo();
}


