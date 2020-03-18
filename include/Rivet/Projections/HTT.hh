// -*- C++ -*-
#ifndef RIVET_HTT_HH
#define RIVET_HTT_HH

#include "Rivet/Jet.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "../../HEPTopTagger/HEPTopTagger.hh"
//#include <functional>

using namespace HEPTopTagger;

namespace Rivet {


  class HTT : public Projection
  {
    private:

        HEPTopTagger::HEPTopTagger _tagger;

        std::vector<double> _topjets;
        bool _do_qjets;

        double _mass_drop_treshold;
        double _max_subjet_mass;

        unsigned _filtering_n;
        double _filtering_R;
        double _filtering_minpT_subjet;
//        JetAlgorithm _filtering_jetalg;
//
//        JetAlgorithm _reclustering_jetalg;


    public:

        HTT(const JetAlg& jetalg);
        DEFAULT_RIVET_PROJ_CLONE(HTT);

        HEPTopTagger::HEPTopTagger* GetNewTagger()
        {
          return &(_tagger);
        }

        HEPTopTagger::HEPTopTagger tagger() {return _tagger;}
        
        void init_jet(const fastjet::PseudoJet& jet) {_tagger(jet);}

        void calc(const Jets& jets);
        
        void Reset();
        
    protected:
    
        void project(const Event& e);
        
        /// Compare projections.
        CmpState compare(const Projection& p) const;

  };


}

#endif
