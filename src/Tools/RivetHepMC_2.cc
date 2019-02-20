// -*- C++ -*-

#include "Rivet/Tools/RivetHepMC.hh"

namespace Rivet{
  
  const Relatives Relatives::PARENTS     = HepMC::parents;
  const Relatives Relatives::CHILDREN    = HepMC::children;
  const Relatives Relatives::ANCESTORS   = HepMC::ancestors;
  const Relatives Relatives::DESCENDANTS = HepMC::descendants;
  
  namespace HepMCUtils{
    
    std::vector<ConstGenParticlePtr> particles(ConstGenEventPtr ge){
      std::vector<ConstGenParticlePtr> result;
      for(GenEvent::particle_const_iterator pi = ge->particles_begin(); pi != ge->particles_end(); ++pi){
        result.push_back(*pi);
      }
      return result;
    }
    
    std::vector<ConstGenParticlePtr> particles(const GenEvent *ge){
      std::vector<ConstGenParticlePtr> result;
      for(GenEvent::particle_const_iterator pi = ge->particles_begin(); pi != ge->particles_end(); ++pi){
        result.push_back(*pi);
      }
      return result;
    }
   
    std::vector<ConstGenVertexPtr> vertices(ConstGenEventPtr ge){
      std::vector<ConstGenVertexPtr> result;
      for(GenEvent::vertex_const_iterator vi = ge->vertices_begin(); vi != ge->vertices_end(); ++vi){
        result.push_back(*vi);
      }
      return result;
    }
    
    std::vector<ConstGenVertexPtr> vertices(const GenEvent *ge){
      std::vector<ConstGenVertexPtr> result;
      for(GenEvent::vertex_const_iterator vi = ge->vertices_begin(); vi != ge->vertices_end(); ++vi){
        result.push_back(*vi);
      }
      return result;
    }
    
    std::vector<ConstGenParticlePtr> particles(ConstGenVertexPtr gv, const Relatives &relo){
      std::vector<ConstGenParticlePtr> result;
      /// @todo A particle_const_iterator on GenVertex would be nice...
      // Before HepMC 2.7.0 there were no GV::particles_const_iterators and constness consistency was all screwed up :-/
#if HEPMC_VERSION_CODE >= 2007000
      for (HepMC::GenVertex::particle_iterator pi = gv->particles_begin(relo); pi != gv->particles_end(relo); ++pi)
      result.push_back(*pi);
#else
      HepMC::GenVertex* gv2 = const_cast<HepMC::GenVertex*>(gv);
      for (HepMC::GenVertex::particle_iterator pi = gv2->particles_begin(relo); pi != gv2->particles_end(relo); ++pi)
      result.push_back(const_cast<ConstGenParticlePtr>(*pi));
#endif
      return result;
    }

    std::vector<ConstGenParticlePtr> particles(ConstGenParticlePtr gp, const Relatives &relo){
      ConstGenVertexPtr vtx;
      
      switch(relo){
        case HepMC::parents:
        
        case HepMC::ancestors:
          vtx = gp->production_vertex();
        break;
        
        case HepMC::children:
        
        case HepMC::descendants:
          vtx = gp->end_vertex();
        break;
        
        default:
        
        throw std::runtime_error("Not implemented!");
        break;
      }
      
      return particles(vtx, relo);
    }

    
    
    int uniqueId(ConstGenParticlePtr gp){
      return gp->barcode();
    }
    
    int particles_size(ConstGenEventPtr ge){
      return ge->particles_size();
    }

    int particles_size(const GenEvent *ge){
      return ge->particles_size();
    }
    
    std::vector<ConstGenParticlePtr> beams(const GenEvent *ge){
      pair<HepMC::GenParticle*, HepMC::GenParticle*> beams = ge->beam_particles();
      return std::vector<ConstGenParticlePtr>{beams.first, beams.second};
    }
    
    std::shared_ptr<HepMC::IO_GenEvent> makeReader(std::istream &istr){
      return make_shared<HepMC::IO_GenEvent>(istr);
    }
   
    bool readEvent(std::shared_ptr<HepMC::IO_GenEvent> io, std::shared_ptr<GenEvent> evt){
      if(io->rdstate() != 0) return false;
      if(!io->fill_next_event(evt.get())) return false;
      return true;
    }
    
  }
}
