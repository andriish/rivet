// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PVertex class.
//

#include "Rivet/Projections/PVertex.hh"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"
#include <stdexcept>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Beam.tcc"
#endif

using namespace Rivet;
using std::runtime_error;

void PVertex::project(const Event& e) {
  //vector<GenParticle*> inc = e.genEvent().signal_process_vertex()->listParents();
  cout << "Entered PVertex.cc" << endl;
  const Event* tempEvent = &e; 
  cout << "tempEvent assigned" << endl;
  //thePVertex = *e.genEvent().signal_process_vertex();
  GenEvent tempGenEvent = tempEvent->genEvent();
  cout << "tempGenEvent assigned" << endl;

  cout << "tempGenEvent.signal_process_id()=" << tempGenEvent.signal_process_id() << endl;
  cout << "tempGenEvent.event_number()=" << tempGenEvent.event_number() << endl;
  cout << "tempGenEvent.event_scale()=" << tempGenEvent.event_scale() << endl;
  cout << "tempGenEvent.pdf_info()=" << tempGenEvent.pdf_info() << endl;
  cout << "tempGenEvent.alphaQCD()=" << tempGenEvent.alphaQCD() << endl;
  cout << "tempGenEvent.alphaQED()=" << tempGenEvent.alphaQED() << endl;


  GenEvent::particle_const_iterator pp; // = ge.particles_begin();
  for (pp = tempGenEvent.particles_begin(); pp!=tempGenEvent.particles_end(); ++pp) {
    cout << "Particle ID=" << (*pp)->pdg_id() << " status=" << (*pp)->status(); //<< endl;
    GenVertex* prodVertex = (*pp)->production_vertex();
    //cout << "Vertex z=" << prodVertex->position().z() << endl;
  }



  for (GenEvent::vertex_const_iterator vit = tempGenEvent.vertices_begin(); vit!= tempGenEvent.vertices_end(); ++vit) {
    //pp = (*vit)->particles_begin();
    //GenEvent::particle_iterator pit = (*vit)->particle_iterator();
    cout << "vertex.id()=" << (*vit)->id()
      //<< " particle_ID=" << (*pp)->pdg_id()
      //<< " particle_ID=" << (*vit)->particle_iterator()->pdg_id()
	 << " particle_ID=" << (*((*vit)->particles_in_const_begin()))->pdg_id()
	 << " particle_status=" << (*((*vit)->particles_in_const_begin()))->status()
      //<< " production_vertex-barcode=" << (*((*vit)->particles_in_const_begin()))->production_vertex()->barcode()
	   << "  position().px()=" << (*vit)->position().px() << "  py=" << (*vit)->position().py()<< "  pz=" << (*vit)->position().pz()
	   << "  barcode()=" << (*vit)->barcode() << endl;
  }




  thePVertex = *(tempGenEvent.signal_process_vertex());
  cout << "Did: thePVertex = *e.genEvent().signal_process_vertex();" << endl;
  if ( thePVertex.particles_in_size() != 2 ) {
    throw std::runtime_error("Wrong number of Primary Vertex particles.");
    /// @todo Maybe we should have our own exception classes?
  }

  /*
  cout << thePVertex->position().t() << " ";
  cout << thePVertex->position().x() << " ";
  cout << thePVertex->position().y() << " ";
  cout << thePVertex->position().z() << endl;
  */



}

int PVertex::compare(const Projection &) const {
  return 0;
}
