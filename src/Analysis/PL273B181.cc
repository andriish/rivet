// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/PL273B181.hh"
#include "HepPDT/ParticleID.hh"

#include "AIDA/IAxis.h"
#include "AIDA/IHistogram1D.h"

using namespace AIDA;
using namespace Rivet;
using namespace HepMC;
using namespace std;


void PL273B181::init() {
  // Book histograms
  histChTot_       = bookHistogram1D("TotalChMult","Total charged multiplicity", 25, 1.0, 51.0);
  histSphericity_  = bookHistogram1D("Sphericity", "Event Shape: Sphericity", 8, 0.0, 0.70);
  histAplanarity_  = bookHistogram1D("Aplanarity", "Event Shape: APlanarity", 10, 0.0, 0.09);
  histPlanarity_   = bookHistogram1D("Planarity",  "Event Shape: Planarity", 16, 0.0, 0.70);
}


// Do the analysis
void PL273B181::analyze(const Event & event) {
  Log& log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;

  // Analyse and print some info
  const Multiplicity& m = event.applyProjection(mult);
  log << Log::INFO << "Total charged multiplicity    = " << m.totalChargedMultiplicity() << endl;

  //Analyse the event shape info
  const Sphericity& s = event.applyProjection(spher);
  log << Log::INFO << "Sphericity    = " << s.sphericity() << endl;
  log << Log::INFO << "Aplanarity    = " << s.aplanarity() << endl;
  log << Log::INFO << "Planarity     = " << s.planarity() << endl;

  // Fill histograms here, and scale them later
  histChTot_->fill(m.totalChargedMultiplicity(), 1.0);
  histSphericity_->fill(s.sphericity(), 1.0);
  histPlanarity_->fill(s.planarity(), 1.0);
  histAplanarity_->fill(s.aplanarity(), 1.0);
  
  // Finished...
  log << Log::DEBUG << "Finished analyzing" << endl;
}

// Finalize
void PL273B181::finalize() { 
  Log& log = getLog();

  double area = 0;
  for (int i=0; i < histSphericity_->axis().bins(); ++i) {
    area += histSphericity_->binHeight(i) * histSphericity_->axis().binWidth(i); 
  }
  log << Log::INFO << "Area under histogram: " << area << endl;


  // Normalize the histogram areas to 1
  if (histChTot_->sumBinHeights()!=0){
    histChTot_->scale(1/histChTot_->sumBinHeights() );
  }
  //histSphericity_->scale(1/histSphericity_->sumBinHeights() );
  histSphericity_->scale(1/area);
  histPlanarity_->scale(1/histPlanarity_->sumBinHeights() ); 
  histAplanarity_->scale(1/histAplanarity_->sumBinHeights() );

  log << Log::INFO << "Sum of sph bin heights after normalization: " 
      << histSphericity_->sumBinHeights() << endl;
}


RivetInfo PL273B181::getInfo() const {
  return Analysis::getInfo() + fsproj.getInfo() +
    mult.getInfo() + spher.getInfo();
}
