// -*- C++ -*-
#ifndef RIVET_HTT_HH
#define RIVET_HTT_HH

#include "Rivet/Jet.hh"
//#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "../../HEPTopTagger/HEPTopTagger.hh"
//#include <functional>

namespace Rivet {


  class HTT : public Projection
  {
    private:

        std::vector<double> _topjets;

    public:

        HTT(const JetAlg& jetalg);
        DEFAULT_RIVET_PROJ_CLONE(HTT);

        void calc(const Jets& jets);
        
        void Reset();
        
    protected:
    
        void project(const Event& e);
        
        /// Compare projections.
        CmpState compare(const Projection& p) const;

  };


}

#endif
