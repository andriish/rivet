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

using namespace HEPTopTagger;
using namespace fastjet;


  class HTT : public Projection
  {
    private:

        std::vector<double> _topjets;

        bool _do_qjets, _do_optimalR, _debug;

        double _mass_drop_treshold, _max_subjet_mass;

        HEPTopTagger::Mode _mode; 
        // (0 = EARLY_MASSRATIO_SORT_MASS)
        // (1 = LATE_MASSRATIO_SORT_MASS)
        // (2 = EARLY_MASSRATIO_SORT_MODDJADE)
        // (3 = LATE_MASSRATIO_SORT_MODDJADE)
        // (4 = TWO_STEP_FILTER)

        double _mtmass, _mwmass, _mtmin, _mtmax, _rmin, _rmax, _fw; 

        double _m23cut, _m13cutmin, _m13cutmax, _minpt_tag;

        unsigned _filtering_n;
        double _filtering_R, _filtering_minpT_subjet;

        JetAlgorithm _filtering_jetalg, _reclustering_jetalg;

        double _zcut;
        double _rcut_factor;

        double _max_fatjet_R; // max should be same with clustering def
        double _min_fatjet_R;
        double _step_R;
        double _optimalR_threshold;

        double _R_filt_optimalR_calc;
        double _N_filt_optimalR_calc;
        double _r_min_exp_function; // =_R_filt_optimalR_calc

        double _optimalR_mmin;
        double _optimalR_mmax;
        double _optimalR_fw;
        double _R_opt_diff;

        double _R_filt_optimalR_pass;
        double _N_filt_optimalR_pass, _R_filt_optimalR_fail, _N_filt_optimalR_fail;

        double _q_zcut, _q_dcut_fctr, _q_exp_min, _q_exp_max, _q_rigidity, _q_truncation_fctr;


    public:

        HTT(const JetAlg& jetalg, 
            unsigned mode=4) : _do_optimalR(1), _do_qjets(0),
                                   _mass_drop_treshold(0.8), _max_subjet_mass(30.),
                                   _mtmass(172.3), _mwmass(80.4), _fw(0.15),
                                   _mtmin(150.), _mtmax(200.), _rmin(0.85*80.4/172.3), _rmax(1.15*80.4/172.3),
                                   _m23cut(0.35), _m13cutmin(0.2), _m13cutmax(1.3), _minpt_tag(200.),
                                   _filtering_n(5), _filtering_R(0.3), _filtering_jetalg(fastjet::cambridge_algorithm), 
                                   _filtering_minpT_subjet(0.), _reclustering_jetalg(fastjet::cambridge_algorithm),
                                   _zcut(0.1), _rcut_factor(0.5),_max_fatjet_R(1.5), 
                                   _min_fatjet_R(0.5), _step_R(0.1), _optimalR_threshold(0.2),
                                   _R_filt_optimalR_calc(0.2), _N_filt_optimalR_calc(10.),
                                   _r_min_exp_function(&_R_filt_optimalR_calc),
                                   _optimalR_mmin(150.), _optimalR_mmax(200.), 
                                   _optimalR_fw(0.175), _R_opt_diff(0.3),
                                   _R_filt_optimalR_pass(0.2), _N_filt_optimalR_pass(5), 
                                   _R_filt_optimalR_fail(0.3), _N_filt_optimalR_fail(3),
                                   _q_zcut(0.1), _q_dcut_fctr(0.5), _q_exp_min(0.), 
                                   _q_exp_max(0.), _q_rigidity(0.1), _q_truncation_fctr(0.); 
        {
            setName("HEPTopTagger");
            declare(jetalg, "Jets");
            set_mode(mode);
        }

        DEFAULT_RIVET_PROJ_CLONE(HTT);

        // Set parameters
        void set_do_qjets(bool m) {_do_qjets=m;}
        void set_do_optimalR(bool m) {_do_optimalR=m;}
        void set_debug(bool m) {_debug=m;}

        void set_mass_drop_treshold(double v) {_mass_drop_treshold=v;}
        void set_max_subjet_mass(double v) {_max_subjet_mass=v;}

        void set_mode(int m) 
        {
            if (m==0)      _mode=HEPTopTagger::EARLY_MASSRATIO_SORT_MASS;
            else if (m==1) _mode=HEPTopTagger::LATE_MASSRATIO_SORT_MASS;
            else if (m==2) _mode=HEPTopTagger::EARLY_MASSRATIO_SORT_MODDJADE;
            else if (m==3) _mode=HEPTopTagger::LATE_MASSRATIO_SORT_MODDJADE;
            else if (m==4) _mode=HEPTopTagger::TWO_STEP_FILTER;
            else 
            {
                MSG_ERROR("Mode has to be between 0 and 4");
                exit(1);
            }
        }

        void set_mtmass(double v) {_mtmass=v;}
        void set_mwmass(double v) {_mwmass=v;}
        void set_mtmin(double v) {_mtmin=v;}
        void set_mtmax(double v) {_mtmax=v;}
        void set_rmin(double v) {_rmin=v;}
        void set_rmax(double v) {_rmax=v;}
        void set_fw(double fw) {_fw=fw;}

        void set_m23cut(double v) {_m23cut=v;}
        void set_m13cutmin(double v) {_m13cutmin=v;}
        void set_m13cutmax(double v) {_m13cutmin=v;}
        void set_minpt_tag(double v) {_minpt_tag=v;}

        void set_filtering_n(unsigned i) {_filtering_n=i;}
        void set_filtering_R(double v) {_filtering_R=v;}
        void set_filtering_minpT_subjet(double v) {_filtering_minpT_subjet=v;}

        void set_filtering_jetalg(JetAlgorithm j) {_filtering_jetalg=j;}
        void set_reclustering_jetalg(JetAlgorithm j) {_reclustering_jetalg=j;}



        void calc(const Jets& jets);
        
        void Reset();
        
    protected:
    
        void project(const Event& e);
        
        /// Compare projections.
        CmpState compare(const Projection& p) const;

  };


}

#endif
