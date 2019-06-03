// -*- C++ -*-
#ifndef RIVET_RivetHepMC_HH
#define RIVET_RivetHepMC_HH

#ifdef ENABLE_HEPMC_3
#include "HepMC3/HepMC3.h"
#include "HepMC3/Relatives.h"
#include "HepMC3/Reader.h"

namespace Rivet{
  namespace RivetHepMC = HepMC3;
  using RivetHepMC::ConstGenParticlePtr;
  using RivetHepMC::ConstGenVertexPtr;
  using RivetHepMC::Relatives;
  using RivetHepMC::ConstGenHeavyIonPtr;
  
  using HepMC_IO_type = RivetHepMC::Reader;

  using PdfInfo = RivetHepMC::GenPdfInfo;
}

#else
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/HeavyIon.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Version.h"
#include "HepMC/GenRanges.h"
#include "HepMC/IO_GenEvent.h"

namespace Rivet{
  namespace RivetHepMC = HepMC;

  
  // HepMC 2.07 provides its own #defines
  typedef const HepMC::GenParticle* ConstGenParticlePtr;
  typedef const HepMC::GenVertex* ConstGenVertexPtr;
  typedef const HepMC::HeavyIon* ConstGenHeavyIonPtr;
  
  /// @brief Replicated the HepMC3 Relatives syntax using HepMC2 IteratorRanges
  /// This is necessary mainly because of capitalisation differences
  class Relatives{
    
    public:
    
    constexpr Relatives(HepMC::IteratorRange relo): _internal(relo){}
    
    constexpr HepMC::IteratorRange operator()() const {return _internal;}
    operator HepMC::IteratorRange() const {return _internal;}
    
    const static Relatives PARENTS;
    const static Relatives CHILDREN;
    const static Relatives ANCESTORS;
    const static Relatives DESCENDANTS;
    
    private:
    const HepMC::IteratorRange _internal;
    
  };
  
  using HepMC_IO_type = HepMC::IO_GenEvent;
  using PdfInfo = RivetHepMC::PdfInfo;

}
  
#endif

#include "Rivet/Tools/RivetSTL.hh"
#include "Rivet/Tools/Exceptions.hh"


namespace Rivet {

  using RivetHepMC::GenEvent;
  using ConstGenEventPtr = std::shared_ptr<const GenEvent>;
  /// @todo Use mcutils?

  namespace HepMCUtils{

    ConstGenParticlePtr              getParticlePtr(const RivetHepMC::GenParticle & gp);
    std::vector<ConstGenParticlePtr> particles(ConstGenEventPtr ge);
    std::vector<ConstGenParticlePtr> particles(const GenEvent *ge);
    std::vector<ConstGenVertexPtr>   vertices(ConstGenEventPtr ge);
    std::vector<ConstGenVertexPtr>   vertices(const GenEvent *ge);
    std::vector<ConstGenParticlePtr> particles(ConstGenVertexPtr gv, const Relatives &relo);
    std::vector<ConstGenParticlePtr> particles(ConstGenParticlePtr gp, const Relatives &relo);
    int uniqueId(ConstGenParticlePtr gp);
    int particles_size(ConstGenEventPtr ge);
    int particles_size(const GenEvent *ge);
    std::pair<ConstGenParticlePtr,ConstGenParticlePtr> beams(const GenEvent *ge);
    std::shared_ptr<HepMC_IO_type> makeReader(std::istream &istr,
                                              std::string * errm = 0);
    bool readEvent(std::shared_ptr<HepMC_IO_type> io,
                   std::shared_ptr<GenEvent> evt);
    void strip(GenEvent & ge,
               const set<long> & stripid = {1, -1, 2, -2, 3,-3, 21});
  }

  inline GenVertexIterRange particles_out(GenVertex* gv) {
    return GenVertexIterRange(gv->particles_out_begin(), gv->particles_out_end());
  }

  #endif


  //////////////////////////


  /// Get the direct parents or all-ancestors of GenParticle @a gp
  inline std::vector<const GenParticle*> particles_in(const GenParticle* gp, HepMC::IteratorRange range=HepMC::ancestors) {
    if (range != HepMC::parents && range != HepMC::ancestors)
      throw UserError("Requested particles_in(GenParticle*) with a non-'in' iterator range");
    if (!gp->production_vertex()) return std::vector<const GenParticle*>();
    #if HEPMC_VERSION_CODE >= 2007000
    return particles(gp->production_vertex(), range);
    #else
    // Before HepMC 2.7.0 the constness consistency of methods and their return types was all screwed up :-/
    std::vector<const GenParticle*> rtn;
    foreach (GenParticle* gp2, particles(gp->production_vertex(), range))
      rtn.push_back( const_cast<const GenParticle*>(gp2) );
    return rtn;
    #endif
  }

  /// Get the direct parents or all-ancestors of GenParticle @a gp
  inline std::vector<GenParticle*> particles_in(GenParticle* gp, HepMC::IteratorRange range=HepMC::ancestors) {
    if (range != HepMC::parents && range != HepMC::ancestors)
      throw UserError("Requested particles_in(GenParticle*) with a non-'in' iterator range");
    return (gp->production_vertex()) ? particles(gp->production_vertex(), range) : std::vector<GenParticle*>();
  }


  /// Get the direct children or all-descendents of GenParticle @a gp
  inline std::vector<const GenParticle*> particles_out(const GenParticle* gp, HepMC::IteratorRange range=HepMC::descendants) {
    if (range != HepMC::children && range != HepMC::descendants)
      throw UserError("Requested particles_out(GenParticle*) with a non-'out' iterator range");
    if (!gp->end_vertex()) return std::vector<const GenParticle*>();
    #if HEPMC_VERSION_CODE >= 2007000
    return particles(gp->end_vertex(), range);
    #else
    // Before HepMC 2.7.0 the constness consistency of methods and their return types was all screwed up :-/
    std::vector<const GenParticle*> rtn;
    foreach (GenParticle* gp2, particles(gp->end_vertex(), range))
      rtn.push_back( const_cast<const GenParticle*>(gp2) );
    return rtn;
    #endif
  }

  /// Get the direct children or all-descendents of GenParticle @a gp
  inline std::vector<GenParticle*> particles_out(GenParticle* gp, HepMC::IteratorRange range=HepMC::descendants) {
    if (range != HepMC::children && range != HepMC::descendants)
      throw UserError("Requested particles_out(GenParticle*) with a non-'out' iterator range");
    return (gp->end_vertex()) ? particles(gp->end_vertex(), range) : std::vector<GenParticle*>();
  }


  /// Get any relatives of GenParticle @a gp
  inline std::vector<const GenParticle*> particles(const GenParticle* gp, HepMC::IteratorRange range=HepMC::ancestors) {
    if (range == HepMC::parents || range == HepMC::ancestors)
      return particles_in(gp, range);
    if (range == HepMC::children || range == HepMC::descendants)
      return particles_in(gp, range);
    throw UserError("Requested particles(GenParticle*) with an unsupported iterator range");
  }


}

#endif
