// -*- C++ -*-

//#include <regex>
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include "Rivet/Tools/Logging.hh"
#include "../Core/zstr/zstr.hpp"

/*namespace {

  inline std::vector<std::string> split(const std::string& input, const std::string& regex) {
    // passing -1 as the submatch index parameter performs splitting
    std::regex re(regex);
    std::sregex_token_iterator
      first{input.begin(), input.end(), re, -1},
      last;
      return {first, last};
  }

}*/

namespace Rivet{
  
  const Relatives Relatives::PARENTS     = HepMC::parents;
  const Relatives Relatives::CHILDREN    = HepMC::children;
  const Relatives Relatives::ANCESTORS   = HepMC::ancestors;
  const Relatives Relatives::DESCENDANTS = HepMC::descendants;
  
  namespace HepMCUtils{
    
    ConstGenParticlePtr getParticlePtr(const RivetHepMC::GenParticle & gp) {
      return &gp;
    }
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
    
    std::pair<ConstGenParticlePtr,ConstGenParticlePtr> beams(const GenEvent *ge){
      return ge->beam_particles();
    }

    std::shared_ptr<HepMC::IO_GenEvent> makeReader(std::string filename,
                                                   std::shared_ptr<std::istream> & istrp,
                                                   std::string *) {
#ifdef HAVE_LIBZ
      if ( filename == "-" )
        istrp = make_shared<Rivet::zstr::istream>(std::cin);
      else
        istrp = make_shared<Rivet::zstr::ifstream>(filename.c_str());
      std::istream & istr = *istrp;
#else
      if ( filename != "-" ) istrp = make_shared<std::ifstream>(filename.c_str());
      std::istream & istr = filename == "-"? std::cin: *istrp;
#endif

      return make_shared<HepMC::IO_GenEvent>(istr);
    }
   
    bool readEvent(std::shared_ptr<HepMC::IO_GenEvent> io, std::shared_ptr<GenEvent> evt){
      if(io->rdstate() != 0) return false;
      if(!io->fill_next_event(evt.get())) return false;
      return true;
    }

    // This functions could be filled with code doing the same stuff as
    // in the HepMC3 version of This file.
    void strip(GenEvent &, const set<long> &) {}

    vector<string> weightNames(const GenEvent & ge) {
      /// reroute the print output to a std::stringstream and process
      /// The iteration is done over a map in hepmc2 so this is safe
      vector<string> ret;

      /// Obtaining weight names using regex probably neater, but regex
      /// is not defined in GCC4.8, which is currently used by Lxplus.
      /// Attempt an alternative solution based on stringstreams:
      std::stringstream stream;
      ge.weights().print(stream);
      std::string pair; // placeholder for subtsring matches
      while (std::getline(stream, pair, ' ')) {
        if ( pair.size() < 2 ) continue;
        pair.erase(pair.begin()); // removes the "(" on the LHS
        pair.pop_back();          // removes the ")" on the RHS
        if (pair.empty())  continue;
        std::stringstream spair(pair);
        vector<string> temp;
        while (std::getline(spair, pair, ',')) {
          temp.push_back(std::move(pair));
        }
        if (temp.size() == 2)  ret.push_back(temp[0]);
      }
      /// Possible future solution based on regex
      /*std::ostringstream stream;
      ge.weights().print(stream);  // Super lame, I know
      string str =  stream.str();

      std::regex re("(([^()]+))"); // Regex for stuff enclosed by parentheses ()
      for (std::sregex_iterator i = std::sregex_iterator(str.begin(), str.end(), re);
           i != std::sregex_iterator(); ++i ) {
        std::smatch m = *i;
        vector<string> temp = ::split(m.str(), "[,]");
        if (temp.size() ==2) {
          // store the default weight based on weight names
          if (temp[0] == "Weight" || temp[0] == "0" || temp[0] == "Default") {
            ret.push_back("");
          } else
            ret.push_back(temp[0]);
        }
      }*/
      return ret;
    }

    pair<double,double> crossSection(const GenEvent & ge) {
      return make_pair(ge.cross_section()->cross_section(),
                       ge.cross_section()->cross_section_error());
    }

    std::valarray<double> weights(const GenEvent & ge) {
      const size_t W = ge.weights().size();
      std::valarray<double> wts(W);
      for (unsigned int iw = 0; iw < W; ++iw)
        wts[iw] = ge.weights()[iw];
      return wts;
    }
  }
}
