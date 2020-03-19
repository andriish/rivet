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

namespace Rivet {

  class HTT : public Projection
  {
    private:

        std::vector<double> _topjets;

        bool _do_qjets=0;
        bool _do_optimalR=1;
        bool _debug=1;

        double _mass_drop_treshold =0.8;
        double _max_subjet_mass=30.;

        HEPTopTagger::Mode _mode=HEPTopTagger::TWO_STEP_FILTER; 
        // (0 = EARLY_MASSRATIO_SORT_MASS)
        // (1 = LATE_MASSRATIO_SORT_MASS)
        // (2 = EARLY_MASSRATIO_SORT_MODDJADE)
        // (3 = LATE_MASSRATIO_SORT_MODDJADE)
        // (4 = TWO_STEP_FILTER)

        double _mtmass = 172.9;
        double _mwmass = 80.379;
        double _mtmin  = 150.;
        double _mtmax  = 200.;
        double _rmin   = 0.85*_mwmass/_mtmass;
        double _rmax   = 1.15*_mwmass/_mtmass;
        double _fw     = 0.15; 

        double _m23cut    = 0.35;
        double _m13cutmin = 0.2;
        double _m13cutmax = 1.3;
        double _minpt_tag = 200.;

        unsigned _filtering_n = 5;
        double _filtering_R   = 0.3;
        double _filtering_minpT_subjet = 0.;

        fastjet::JetAlgorithm _filtering_jetalg    = fastjet::cambridge_algorithm;
        fastjet::JetAlgorithm _reclustering_jetalg = fastjet::cambridge_algorithm;

        double _zcut        = 0.1;
        double _rcut_factor = 0.5;

        double _max_fatjet_R       = 1.5; // max should be same with clustering def
        double _min_fatjet_R       = 0.5;
        double _step_R             = 0.1;
        double _optimalR_threshold = 0.2;

        double _R_filt_optimalR_calc = 0.2;
        double _N_filt_optimalR_calc = 10.;
        double _r_min_exp_function   = _R_filt_optimalR_calc; // =_R_filt_optimalR_calc

        double _optimalR_mmin = 150.;
        double _optimalR_mmax = 200.;
        double _optimalR_fw   = 0.175;
        double _R_opt_diff    = 0.3;

        double _R_filt_optimalR_pass = .2;
        double _N_filt_optimalR_pass = 5.;
        double _R_filt_optimalR_fail = .3;
        double _N_filt_optimalR_fail = 3.;

        double _q_zcut =0.1;
        double _q_dcut_fctr =0.5;
        double _q_exp_min =0.;
        double _q_exp_max = .0;
        double _q_rigidity =.1;
        double _q_truncation_fctr = 0.;


    public:

        HTT(const JetAlg& jetalg, const std::map<std::string,std::string>& options)
        {
            setName("HEPTopTagger");
            declare(jetalg, "Jets");
            Set_Parameters(options);
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

        void set_filtering_jetalg(fastjet::JetAlgorithm j) {_filtering_jetalg=j;}
        void set_reclustering_jetalg(fastjet::JetAlgorithm j) {_reclustering_jetalg=j;}

        void Set_Parameters(const std::map<std::string, std::string>& options);

        void calc(const Jets& jets);
        
        void Reset();
        
    protected:
    
        void project(const Event& e);
        
        /// Compare projections.
        CmpState compare(const Projection& p) const;

  };


}

#endif
