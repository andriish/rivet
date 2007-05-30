// -*- C++ -*-

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
  //cout << "Entered PVertex.cc" << endl;

  //for some reason I can't access the signal_process_vertex, as if Event to GenEevent
  //conversion would not propagate this information (get segmentation violation)
  //so I take just the first vertex in the HepMC list which should be located at the same 
  //z() pisition as the primary vertex
  //@todo needs further investigation

  /*
  const Event* tempEvent = &e; 
  cout << "tempEvent assigned" << endl;
  //thePVertex = *e.genEvent().signal_process_vertex();
  GenEvent tempGenEvent = tempEvent->genEvent();
  cout << "tempGenEvent assigned" << endl;
  */

  //HepMC::GenEvent genevt = e.genEvent();
  const HepMC::GenEvent & genevt = e.genEvent();

  /*
  cout << "tempGenEvent.signal_process_id()=" << tempGenEvent.signal_process_id() << endl;
  cout << "tempGenEvent.event_number()=" << tempGenEvent.event_number() << endl;
  cout << "tempGenEvent.event_scale()=" << tempGenEvent.event_scale() << endl;
  cout << "tempGenEvent.pdf_info()=" << tempGenEvent.pdf_info() << endl;
  cout << "tempGenEvent.alphaQCD()=" << tempGenEvent.alphaQCD() << endl;
  cout << "tempGenEvent.alphaQED()=" << tempGenEvent.alphaQED() << endl;
  */

  //cout << "genevt assigned!" << endl;
  //HepMC::GenVertex* sigvertex = tempGenEvent.signal_process_vertex();
  /*
  HepMC::GenVertex* sigvertex = genevt.signal_process_vertex();
  cout << "sigvertex assigned!" << endl;
  
  cout << "vtx->t=" << sigvertex->position().t() << " ";
  cout << "vtx->x=" << sigvertex->position().x() << " ";
  cout << "vtx->y=" << sigvertex->position().y() << " ";
  cout << "vtx->z=" << sigvertex->position().z() << endl;
  

  HepMC::GenVertex::particles_in_const_iterator pp = sigvertex->particles_in_const_begin();
  cout << "Signal Vertex Size (In) = " << sigvertex->particles_in_size() 
       << ", (Out) = " << sigvertex->particles_out_size()
       << ", Vertex ID=" << sigvertex->id() << ", barcode=" 
       << sigvertex->barcode() << endl;

  HepMC::GenParticle* first = *pp;
  ++pp;
  HepMC::GenParticle* second = *pp;
  
  cout << "1st particle: " << first->pdg_id() << endl;
  cout << "2nd particle: " << second->pdg_id() << endl;
  */

  /*
  GenEvent::particle_const_iterator pp; // = ge.particles_begin();
  for (pp = tempGenEvent.particles_begin(); pp!=tempGenEvent.particles_end(); ++pp) {
    //cout << "Particle ID=" << (*pp)->pdg_id() << " status=" << (*pp)->status(); //<< endl;
    //GenVertex* prodVertex = (*pp)->production_vertex();
    //cout << "Vertex z=" << prodVertex->position().z() << endl;
  }
  */

  /*
  //for (GenEvent::vertex_const_iterator vit = tempGenEvent.vertices_begin(); vit!= tempGenEvent.vertices_end(); ++vit) {
  for (GenEvent::vertex_const_iterator vit = genevt.vertices_begin(); vit!= genevt.vertices_end(); ++vit) {
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
  */


  //cout << "copy the PV!" << endl;
  //thePVertex = *(tempGenEvent.signal_process_vertex());
  //thePVertex = *(genevt.signal_process_vertex());
  GenEvent::vertex_const_iterator vit = genevt.vertices_begin();
  /*
  for (vit = genevt.vertices_begin(); 
       GenVertex(**vit).particles_in_size()!=2 
	 && vit != genevt.vertices_end(); 
       ++vit) {
    cout << "vit: size=" << GenVertex(**vit).particles_in_size() << endl;
    }
  */

  thePVertex = GenVertex(**vit);



  //cout << "Did: thePVertex = *e.genEvent().signal_process_vertex();" << endl;
  //cout << "Did: thePVertex = GenVertex(**vit) size()=" << thePVertex.particles_in_size() << endl;

  //if ( thePVertex.particles_in_size() != 2 ) {
  if ( thePVertex.particles_in_size() == 0 ) {
    throw std::runtime_error("Wrong number of Primary Vertex particles.");
    /// @todo Maybe we should have our own exception classes?
  }

  /* 
  cout << "PVertex t=" << thePVertex.position().t() << " ";
  cout << "PVertex x=" << thePVertex.position().x() << " ";
  cout << "PVertex y=" << thePVertex.position().y() << " ";
  cout << "PVertex z=" << thePVertex.position().z() << endl;
  */



}

int PVertex::compare(const Projection &) const {
  return 0;
}
