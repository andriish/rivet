// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/HepEx0701051.hh"

#include "Rivet/RivetAIDA.hh"
using namespace Rivet;

#include "AIDA/IHistogram1D.h"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////


// Initialization definition, in this case simply books a 1D histogram
void HepEx0701051::init() {
  vector<double> binedges; 
  binedges.push_back(54);
  binedges.push_back(62);
  binedges.push_back(72);
  binedges.push_back(83);
  binedges.push_back(96);
  binedges.push_back(110);
  binedges.push_back(127);
  binedges.push_back(146);
  binedges.push_back(169);
  binedges.push_back(195);
  binedges.push_back(224);
  binedges.push_back(259);
  binedges.push_back(298);
  binedges.push_back(344);
  binedges.push_back(396);
  binedges.push_back(457);
  binedges.push_back(527);
  binedges.push_back(700); 
  
  pHistogramObject1 = bookHistogram1D("eta less than 0.1", "eta less than 0.1", binedges);
  pHistogramObject2 = bookHistogram1D("eta between 0.1 and 0.7", "eta between 0.1 and 0.7", binedges);
  pHistogramObject3 = bookHistogram1D("eta between 0.7 and 1.1", "eta between 0.7 and 1.1", binedges);
  pHistogramObject4 = bookHistogram1D("eta between 1.1 and 1.6", "eta between 1.1 and 1.6", binedges);
  pHistogramObject5 = bookHistogram1D("eta between 1.6 and 2.1", "eta between 1.6 and 2.1", binedges);
  pHistogramObject1N = bookHistogram1D("eta less than 0.1 N", "eta less than 0.1 N", binedges);
  pHistogramObject2N = bookHistogram1D("eta between 0.1 and 0.7 N", "eta between 0.1 and 0.7 N", binedges);
  pHistogramObject3N = bookHistogram1D("eta between 0.7 and 1.1 N", "eta between 0.7 and 1.1 N", binedges);
  pHistogramObject4N = bookHistogram1D("eta between 1.1 and 1.6 N", "eta between 1.1 and 1.6 N", binedges);
  pHistogramObject5N = bookHistogram1D("eta between 1.6 and 2.1 N", "eta between 1.6 and 2.1 N", binedges);
}


// Do the analysis
void HepEx0701051::analyze(const Event& event) {
  Log log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;

  double binEdges[18] = {54, 62, 72, 83, 96, 110, 127, 146, 169, 195, 224, 259, 298, 344, 396, 457, 527, 700};
  
  // Declare projections here
  FinalState fsproj;
  KtJets ktproj(fsproj, 4, 2, 1, 0.7);
  
  // Apply projections
  const KtJets& kt = event.applyProjection(ktproj);
  
  // Extract data from projected objects
  vector<KtJet::KtLorentzVector> jetList = kt.getJets();

  // Fill histograms
  for (std::vector<KtJet::KtLorentzVector>::iterator j = jetList.begin(); j != jetList.end(); j++) {
        
    if (j->perp() >= 54) {
      
      // ITERATE OVER EACH BIN
      for (int i=0; i<17; i++) {
        
        // ETA LESS THAN 0.1
        if (sqrt(j->eta()*j->eta()) < 0.1) {
          if ((j->perp() >= binEdges[i]) && (j->perp() < binEdges[i+1])){
            pHistogramObject1->fill(j->perp());
          }
        }
	
	// ETA BETWEEN 0.1 AND 0.7
	
	if (sqrt(j->eta()*j->eta()) > 0.1 && sqrt(j->eta()*j->eta()) < 0.7){
	  if ((j->perp() >= binEdges[i]) && (j->perp() < binEdges[i+1])){
	    pHistogramObject2->fill(j->perp());
	  }
	}
	
	// ETA BETWEEN 0.7 AND 1.1
	
	if (sqrt(j->eta()*j->eta()) > 0.7 && sqrt(j->eta()*j->eta()) < 1.1){
	  if ((j->perp() >= binEdges[i]) && (j->perp() < binEdges[i+1])){
	    pHistogramObject3->fill(j->perp());
	  }
	}
	
	// ETA BETWEEN 1.1 AND 1.6
	
	if (sqrt(j->eta()*j->eta()) > 1.1 && sqrt(j->eta()*j->eta()) < 1.6){
	  if ((j->perp() >= binEdges[i]) && (j->perp() < binEdges[i+1])){
	    pHistogramObject4->fill(j->perp());
	  }
	}
	
	// ETA BETWEEN 1.6  AND 2.1
	
	if (sqrt(j->eta()*j->eta()) > 1.6  && sqrt(j->eta()*j->eta()) < 2.1) {
	  if ((j->perp() >= binEdges[i]) && (j->perp() < binEdges[i+1])){
	    pHistogramObject5->fill(j->perp());
	  }
	}
      }
    }  
  }
  
  // Finished
  log << Log::DEBUG << "Finished analyzing" << endl;
}


// Finalize
void HepEx0701051::finalize() {
  // luminosity is calculated from MC generator
  double luminosity = 1;
  double binValues[17] = {58, 68, 77.5, 89.5, 103, 118.5, 136.5, 157.5, 182, 209.5, 241.5, 278.5, 321, 370, 426.5, 492, 613.5};
  double binWidths[17] = {8, 10, 11, 13, 14, 17, 19, 23, 26, 29, 35, 39, 46, 52, 61, 70, 173};
  int a = 0;
  for (int i = 0; i < 17; ++i) {
    pHistogramObject1N->fill(binValues[i], (1 / (binWidths[i]*0.1*luminosity)) * pHistogramObject1->binEntries(i));
    pHistogramObject2N->fill(binValues[i], (1 / (binWidths[i]*0.6*luminosity)) * pHistogramObject2->binEntries(i));
    pHistogramObject3N->fill(binValues[i], (1 / (binWidths[i]*0.4*luminosity)) * pHistogramObject3->binEntries(i));
    pHistogramObject4N->fill(binValues[i], (1 / (binWidths[i]*0.5*luminosity)) * pHistogramObject4->binEntries(i));
    pHistogramObject5N->fill(binValues[i], (1 / (binWidths[i]*0.5*luminosity)) * pHistogramObject5->binEntries(i));
    a = a + pHistogramObject1->binEntries(i);
  }
}
