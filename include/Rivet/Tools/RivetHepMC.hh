// -*- C++ -*-
#ifndef RIVET_RivetHepMC_HH
#define RIVET_RivetHepMC_HH

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Version.h"
#include "HepMC/IO_GenEvent.h"
#include "Rivet/Tools/RivetSTL.hh"
#include "Rivet/Tools/Exceptions.hh"

#if HEPMC_VERSION_CODE >= 3000000
#include "HepMC/HepMC.h"
/// @todo Will this be included in the HepMC/HepMC.h header eventually?
#include "HepMC/Relatives.h"
#else
#include "HepMC/GenRanges.h"
#endif

namespace Rivet {


  /// @todo Use mcutils?

  using HepMC::GenEvent;
  /// @todo should we ever refer to a GenParticle other than by ptr?  Who will own the gp?
  //using HepMC::GenParticle;
  /// @todo should we ever refer to a GenVertex other than by ptr?  Who will own the gv?
  //using HepMC::GenVertex;

  #if HEPMC_VERSION_CODE >= 3000000
  //using HepMC::GenParticlePtr;
  using HepMC::ConstGenParticlePtr;
  //using HepMC::GenVertexPtr;
  using HepMC::ConstGenVertexPtr;
  using HepMC::Relatives;
  #elif HEPMC_VERSION_CODE >= 2007000
  // HepMC 2.07 provides its own #defines
  #else
  #define ConstGenParticlePtr const GenParticle*
  #define ConstGenVertexPtr const GenVertex*
  #define Relatives HepMC::IteratorRange
  #endif

  std::vector<ConstGenParticlePtr> particles(const GenEvent *ge);
  std::vector<ConstGenVertexPtr> vertices(const GenEvent *ge);
  std::vector<ConstGenParticlePtr> particles(ConstGenVertexPtr gv, const Relatives &relo);
  std::vector<ConstGenParticlePtr> particles(ConstGenParticlePtr gp, const Relatives &relo);


  #if HEPMC_VERSION_CODE >= 3000000

  /// @name Accessors from GenEvent
  //@{
/*
  inline std::vector<GenParticlePtr> particles(const GenEvent* ge) {
    assert(ge);
    return ge->particles();
    // std::vector<GenParticlePtr> rtn;
    // for (const GenParticlePtr p : ge->particles()) rtn.push_back(p);
    // return rtn;
  }

  inline std::vector<GenVertexPtr> vertices(const GenEvent* ge) {
    assert(ge);
    return ge->vertices();
    // std::vector<GenVertexPtr> rtn;
    // for (const GenVertexPtr v : ge->vertices()) rtn.push_back(v);
    // return rtn;
  }
*/
  //@}


  /// @name Accessors from GenVertex
  //@{

  /// @todo are these really necessary? Why not call GenVertex::particles_[in, out] directly?
  //inline const vector<GenParticlePtr>& particles_in(const GenVertexPtr& gv) { return gv->particles_in(); }
  // inline vector<GenParticlePtr>& particles_in(GenVertexPtr& gv) { return gv->particles_in(); }

  //inline const vector<GenParticlePtr>& particles_out(const GenVertexPtr& gv) { return gv->particles_out(); }
  // inline vector<GenParticlePtr>& particles_out(GenVertexPtr& gv) { return gv->particles_out(); }
/*
  inline std::vector<GenParticlePtr> particles(const GenVertexPtr& gv, HepMC::IteratorRange range) {
    return HepMC::FindParticles(gv, range).results();
  }
*/
  // /// Get the direct parents or all-ancestors of GenParticle @a gp
  // inline std::vector<GenParticlePtr> particles_in(GenParticlePtr gp, HepMC::IteratorRange range=HepMC::ancestors) {
  //   assert(gp);
  //   if (range != HepMC::parents && range != HepMC::ancestors)
  //     throw UserError("Requested particles_in(GenParticlePtr) with a non-'in' iterator range");
  //   return particles(gp->production_vertex(), range);
  // }
  // /// Get the direct children or all-descendents of GenParticle @a gp
  // inline std::vector<GenParticlePtr> particles_out(GenParticlePtr gp, HepMC::IteratorRange range=HepMC::ancestors) {
  //   assert(gp);
  //   if (range != HepMC::children && range != HepMC::descendants)
  //     throw UserError("Requested particles_out(GenParticlePtr) with a non-'out' iterator range");
  //   return particles(gp->production_vertex(), range);
  // }

  //@}


  /// @name Accessors from GenParticle
  //@{

  /// Get any relatives of GenParticle @a gp
  /*
  inline std::vector<GenParticlePtr> particles(GenParticlePtr gp, HepMC::IteratorRange range) {
    return HepMC::FindParticles(gp, range).results();
  }*/
  //@}


  #else


  /// @name Accessors from GenEvent
  //@{

  inline std::vector<const GenParticlePtr> particles(const GenEvent* ge) {
    assert(ge);
    std::vector<const GenParticlePtr> rtn;
    #if HEPMC_VERSION_CODE >= 2007000
    for (const GenParticlePtr p : ge->particles()) rtn.push_back(p);
    #else
    for (GenEvent::particle_const_iterator pi = ge->particles_begin(); pi != ge->particles_end(); ++pi)
      rtn.push_back(*pi);
    #endif
    return rtn;
  }
  inline std::vector<GenParticlePtr> particles(GenEvent* ge) {
    assert(ge);
    std::vector<GenParticlePtr> rtn;
    #if HEPMC_VERSION_CODE >= 2007000
    for (GenParticlePtr p : ge->particles()) rtn.push_back(p);
    #else
    for (GenEvent::particle_iterator pi = ge->particles_begin(); pi != ge->particles_end(); ++pi)
      rtn.push_back(*pi);
    #endif
    return rtn;
  }


  inline std::vector<const GenVertexPtr> vertices(const GenEvent* ge) {
    assert(ge);
    std::vector<const GenVertexPtr> rtn;
    #if HEPMC_VERSION_CODE >= 2007000
    for (const GenVertexPtr v : ge->vertices()) rtn.push_back(v);
    #else
    for (GenEvent::vertex_const_iterator vi = ge->vertices_begin(); vi != ge->vertices_end(); ++vi)
      rtn.push_back(*vi);
    #endif
    return rtn;
  }
  inline std::vector<GenVertexPtr> vertices(GenEvent* ge) {
    assert(ge);
    std::vector<GenVertexPtr> rtn;
    #if HEPMC_VERSION_CODE >= 2007000
    for (const GenVertexPtr v : ge->vertices()) rtn.push_back(v);
    #else
    for (GenEvent::vertex_iterator vi = ge->vertices_begin(); vi != ge->vertices_end(); ++vi)
      rtn.push_back(*vi);
    #endif
    return rtn;
  }

  //@}


  /// @name Accessors from GenVertex
  //@{

  inline std::vector<const GenParticlePtr> particles(const GenVertexPtr gv, HepMC::IteratorRange range=HepMC::relatives) {
    std::vector<const GenParticlePtr> rtn;
    /// @todo A particle_const_iterator on GenVertex would be nice...
    // Before HepMC 2.7.0 there were no GV::particles_const_iterators and constness consistency was all screwed up :-/
    #if HEPMC_VERSION_CODE >= 2007000
    for (const GenParticlePtr p : gv->particles(range)) rtn.push_back(p);
    #else
    GenVertexPtr gv2 = const_cast<GenVertexPtr>(gv);
    for (GenVertex::particle_iterator pi = gv2->particles_begin(range); pi != gv2->particles_end(range); ++pi)
      rtn.push_back(const_cast<const GenParticlePtr>(*pi));
    #endif
    return rtn;
  }
  inline std::vector<GenParticlePtr> particles(GenVertexPtr gv, HepMC::IteratorRange range=HepMC::relatives) {
    std::vector<GenParticlePtr> rtn;
    for (GenVertex::particle_iterator pi = gv->particles_begin(range); pi != gv->particles_end(range); ++pi)
      rtn.push_back(*pi);
    return rtn;
  }


  // Get iterator ranges as wrapped begin/end pairs
  /// @note GenVertex _in and _out iterators are actually, secretly the same types *sigh*
  struct GenVertexIterRangeC {
    typedef vector<GenParticlePtr>::const_iterator genvertex_particles_const_iterator;
    GenVertexIterRangeC(const genvertex_particles_const_iterator& begin, const genvertex_particles_const_iterator& end)
      : _begin(begin), _end(end) {  }
    const genvertex_particles_const_iterator& begin() { return _begin; }
    const genvertex_particles_const_iterator& end() { return _end; }
  private:
    const genvertex_particles_const_iterator _begin, _end;
  };

  inline GenVertexIterRangeC particles_in(const GenVertexPtr gv) {
    return GenVertexIterRangeC(gv->particles_in_const_begin(), gv->particles_in_const_end());
  }

  inline GenVertexIterRangeC particles_out(const GenVertexPtr gv) {
    return GenVertexIterRangeC(gv->particles_out_const_begin(), gv->particles_out_const_end());
  }


  #if HEPMC_VERSION_CODE >= 2007000

  // Get iterator ranges as wrapped begin/end pairs
  /// @note GenVertex _in and _out iterators are actually, secretly the same types *sigh*
  struct GenVertexIterRange {
    typedef vector<GenParticlePtr>::iterator genvertex_particles_iterator;
    GenVertexIterRange(const genvertex_particles_iterator& begin, const genvertex_particles_iterator& end)
      : _begin(begin), _end(end) {  }
    const genvertex_particles_iterator& begin() { return _begin; }
    const genvertex_particles_iterator& end() { return _end; }
  private:
    const genvertex_particles_iterator _begin, _end;
  };

  inline GenVertexIterRange particles_in(GenVertexPtr gv) {
    return GenVertexIterRange(gv->particles_in_begin(), gv->particles_in_end());
  }

  inline GenVertexIterRange particles_out(GenVertexPtr gv) {
    return GenVertexIterRange(gv->particles_out_begin(), gv->particles_out_end());
  }

  #endif


  /// @name Accessors from GenParticle
  //@{

  /// Get the direct parents or all-ancestors of GenParticle @a gp
  inline std::vector<const GenParticlePtr> _particles_in(const GenParticlePtr gp, HepMC::IteratorRange range=HepMC::ancestors) {
    assert(gp);
    if (range != HepMC::parents && range != HepMC::ancestors)
      throw UserError("Requested particles_in(GenParticlePtr) with a non-'in' iterator range");
    if (!gp->production_vertex()) return std::vector<const GenParticlePtr>();
    #if HEPMC_VERSION_CODE >= 2007000
    return particles(gp->production_vertex(), range);
    #else
    // Before HepMC 2.7.0 the constness consistency of methods and their return types was all screwed up :-/
    std::vector<const GenParticlePtr> rtn;
    for (GenParticlePtr gp2 : particles(gp->production_vertex(), range))
      rtn.push_back( const_cast<const GenParticlePtr>(gp2) );
    return rtn;
    #endif
  }
  /// Get the direct parents or all-ancestors of GenParticle @a gp
  inline std::vector<GenParticlePtr> _particles_in(GenParticlePtr gp, HepMC::IteratorRange range=HepMC::ancestors) {
    assert(gp);
    if (range != HepMC::parents && range != HepMC::ancestors)
      throw UserError("Requested particles_in(GenParticlePtr) with a non-'in' iterator range");
    return (gp->production_vertex()) ? particles(gp->production_vertex(), range) : std::vector<GenParticlePtr>();
  }


  /// Get the direct children or all-descendents of GenParticle @a gp
  inline std::vector<const GenParticlePtr> _particles_out(const GenParticlePtr gp, HepMC::IteratorRange range=HepMC::descendants) {
    assert(gp);
    if (range != HepMC::children && range != HepMC::descendants)
      throw UserError("Requested particles_out(GenParticlePtr) with a non-'out' iterator range");
    if (!gp->end_vertex()) return std::vector<const GenParticlePtr>();
    #if HEPMC_VERSION_CODE >= 2007000
    return particles(gp->end_vertex(), range);
    #else
    // Before HepMC 2.7.0 the constness consistency of methods and their return types was all screwed up :-/
    std::vector<const GenParticlePtr> rtn;
    foreach (GenParticlePtr gp2, particles(gp->end_vertex(), range))
      rtn.push_back( const_cast<const GenParticlePtr>(gp2) );
    return rtn;
    #endif
  }
  /// Get the direct children or all-descendents of GenParticle @a gp
  inline std::vector<GenParticlePtr> _particles_out(GenParticlePtr gp, HepMC::IteratorRange range=HepMC::descendants) {
    assert(gp);
    if (range != HepMC::children && range != HepMC::descendants)
      throw UserError("Requested particles_out(GenParticlePtr) with a non-'out' iterator range");
    return (gp->end_vertex()) ? particles(gp->end_vertex(), range) : std::vector<GenParticlePtr>();
  }


  /// Get any relatives of GenParticle @a gp
  inline std::vector<const GenParticlePtr> particles(const GenParticlePtr gp, HepMC::IteratorRange range) {
    if (range == HepMC::parents || range == HepMC::ancestors)
      return _particles_in(gp, range);
    if (range == HepMC::children || range == HepMC::descendants)
      return _particles_out(gp, range);
    throw UserError("Requested particles(const GenParticlePtr) with an unsupported iterator range");
  }
  /// Get any relatives of GenParticle @a gp
  inline std::vector<GenParticlePtr> particles(GenParticlePtr gp, HepMC::IteratorRange range) {
    if (range == HepMC::parents || range == HepMC::ancestors)
      return _particles_in(gp, range);
    if (range == HepMC::children || range == HepMC::descendants)
      return _particles_out(gp, range);
    throw UserError("Requested particles(GenParticlePtr) with an unsupported iterator range");
  }

  #endif


}

#endif
