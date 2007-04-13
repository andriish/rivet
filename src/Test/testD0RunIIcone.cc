
#include <list.h>

#include "energycluster/ILConeAlgorithm.hpp"
#include "HepEntity.h"



int main() {

  
  HepEntity el;
  list<const HepEntity*> *ensemble = new list<const HepEntity*>;


  //fill some test 4-vectors with E, px, py, pz
  el.Fill(100., 25., 25., 25.);
  ensemble->push_back(new HepEntity(el));
  el.Fill(105., 20., 30., 30.);
  ensemble->push_back(new HepEntity(el));
  el.Fill(60., 20., 20., 20.);
  ensemble->push_back(new HepEntity(el));
  el.Fill(95., 65., 10., 20.);
  ensemble->push_back(new HepEntity(el));
  
  el.Fill(110., 25., -25., -25.);
  ensemble->push_back(new HepEntity(el));
  el.Fill(100., 23., -25., -25.);
  ensemble->push_back(new HepEntity(el));
  el.Fill(101., 25., -20., -25.);
  ensemble->push_back(new HepEntity(el));
  el.Fill(102., 25., -25., -23.);
  ensemble->push_back(new HepEntity(el));
  
  //print information about initial 4-vectors
  cout << "list->size()=" << ensemble->size() << endl;
  int i=1;
  for (list<const HepEntity*>::iterator it = ensemble->begin(); it != ensemble->end(); ++it) {
    cout << "4-vector " << i++ << " : E=" << (*it)->E << " pT=" << (*it)->pT() << " y=" << (*it)->y() << " phi=" << (*it)->phi() << endl; 
  }

  //initialize D0RunII cone algorithm
  float cone_radius = 0.5;
  float min_jet_Et = 0.;
  float split_ratio = 0.5;

  float far_ratio=0.5;
  float Et_min_ratio=0.5;
  bool kill_duplicate=true;
  float duplicate_dR=0.005; 
  float duplicate_dPT=0.01; 
  float search_factor=1.0; 
  float pT_min_leading_protojet=0.; 
  float pT_min_second_protojet=0.;
  int merge_max=10000; 
  float pT_min_nomerge=0.;

  ILConeAlgorithm<HepEntity> 
    ilegac(cone_radius, min_jet_Et, split_ratio,
	   far_ratio, Et_min_ratio, kill_duplicate, duplicate_dR, 
	   duplicate_dPT, search_factor, pT_min_leading_protojet, 
	   pT_min_second_protojet, merge_max, pT_min_nomerge);
 
  float Item_ET_Threshold = 0.;
  //float* Item_ET_Threshold_ptr = &Item_ET_Threshold;

  list<HepEntity> jets;
  ilegac.makeClusters(jets, *ensemble, Item_ET_Threshold);

  list<HepEntity>::iterator it;
  cout << "Number of jets = " << jets.size() << endl;
  for (it=jets.begin(); it!=jets.end(); ++it) {
    cout << "jet: E=" << (*it).E << " pT=" << (*it).pT() << " y=" << (*it).y() << " phi=" << (*it).phi() << endl; 
  }

  return 0;

}
