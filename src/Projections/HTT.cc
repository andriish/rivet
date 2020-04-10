// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/HTT.hh"

namespace Rivet {

using namespace std;


//HTT::HTT(const JetAlg& jetalg, const std::map<std::string,std::string>& options)
//{
//    setName("HEPTopTagger");
//    declare(jetalg, "Jets");
//    Set_Parameters(options);
//}

CmpState HTT::compare(const Projection& p) const
{
  /// @todo Needs to depend on params/options and cuts
    return CmpState::EQ;
}

void HTT::Reset()
{
    _topjets.clear();
}

void HTT::calc(const Jets& xjets) {

  /// @todo Remove hard coding!
  Jets jets = select(xjets, Cuts::pT > 200*GeV);
  /// @todo Use range-for syntax: cleaner and more hackable
  for (unsigned i=0; i<jets.size();i++)
    {
        HEPTopTagger::HEPTopTagger tagger(jets[i].pseudojet());

        MSG_INFO("Constituents : " << jets[i].pseudojet().constituents().size());
        // Set parameters
        tagger.set_mass_drop_threshold(_mass_drop_treshold);
        tagger.set_max_subjet_mass(_max_subjet_mass);

        tagger.set_filtering_n(_filtering_n);
        tagger.set_filtering_R(_filtering_R);
        tagger.set_filtering_minpt_subjet(_filtering_minpT_subjet);
        tagger.set_filtering_jetalgorithm(_filtering_jetalg);

        tagger.set_reclustering_jetalgorithm(_reclustering_jetalg);

        // How to select among candidates
        tagger.set_mode(_mode);
        tagger.set_mt(_mtmass);
        tagger.set_mw(_mwmass);
        // Requirements to accept a candidate
        tagger.set_top_mass_range(_mtmin, _mtmax);
        tagger.set_fw(_fw);
        tagger.set_mass_ratio_range(_rmin,_rmax);
        tagger.set_mass_ratio_cut(_m23cut,_m13cutmin,_m13cutmax);
        tagger.set_top_minpt(_minpt_tag);

        tagger.set_optimalR_max(_max_fatjet_R);
        tagger.set_optimalR_min(_min_fatjet_R);
        tagger.set_optimalR_step(_step_R);
        tagger.set_optimalR_threshold(_optimalR_threshold);

        tagger.set_filtering_optimalR_calc_R(_R_filt_optimalR_calc);
        tagger.set_filtering_optimalR_calc_n(_N_filt_optimalR_calc);
//        tagger.set_optimalR_calc_fun(_r_min_exp_function);

        tagger.set_optimalR_type_top_mass_range(_optimalR_mmin, _optimalR_mmax);
        tagger.set_optimalR_type_fw(_optimalR_fw);
        tagger.set_optimalR_type_max_diff(_R_opt_diff);

        tagger.set_filtering_optimalR_pass_R(_R_filt_optimalR_pass);
        tagger.set_filtering_optimalR_pass_n(_N_filt_optimalR_pass);
        tagger.set_filtering_optimalR_fail_R(_R_filt_optimalR_fail);
        tagger.set_filtering_optimalR_fail_n(_N_filt_optimalR_fail);

        tagger.set_pruning_zcut(_zcut);
        tagger.set_pruning_rcut_factor(_rcut_factor);

        tagger.set_debug(_debug);
        tagger.do_qjets(_do_qjets);
        tagger.set_qjets(_q_zcut,_q_dcut_fctr,_q_exp_min,
                         _q_exp_max,_q_rigidity,_q_truncation_fctr);

        // Run the tagger
        tagger.run();
        if (tagger.is_tagged())
            _topjets.push_back(jets[i]);
        if (_debug)
        {
            tagger.get_info();
            tagger.get_setting();
            MSG_INFO("Maybe top: " << tagger.is_maybe_top());
            // Look at output if we have a tag:
            if (tagger.is_tagged()){
                MSG_INFO("Input fatjet: " << i << "  pT = " << jets[i].perp());
                MSG_INFO("Output: pT = " << tagger.t().perp() << " Mass = " << tagger.t().m() << " f_rec = " << tagger.f_rec());
            } else {
                MSG_INFO("Not tagged");
            }
        //MSG_INFO("Jet PT = " << jet.pT()/GeV);
        }
    }
}

void HTT::project(const Event& e) {
    const Jets jets = applyProjection<JetAlg>(e, "Jets").jets();
    Reset();
    calc(jets);
}



void HTT::Set_Parameters(const std::map<std::string, std::string>& options)
{
      // Loop over options
      for (std::map<std::string,std::string>::const_iterator
           it=options.begin();it!=options.end();it++)
        {
            std::string key = it->first;
            std::for_each(key.begin(), key.end(), [](char & c)
            { c = ::tolower(c); });

            if (key == "topmass")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                if (tmp>0.) _mtmass = tmp;
            }
            else if (key == "wmass")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                if (tmp>0.) _mwmass = tmp;
            }
            else if (key == "do_qjet")
            {
                bool tmp=0;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _do_qjets = tmp;
            }
            else if (key == "do_optimalr")
            {
                bool tmp=1;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _do_optimalR = tmp;
            }
            else if (key == "debug")
            {
                bool tmp=0;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _debug = tmp;
            }
            else if (key == "massdroptreshold" || key == "mass_drop_treshold")
            {
                double tmp=0;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _mass_drop_treshold = tmp;
            }
            else if (key == "max_subjet_mass" || key == "maxsubjetmass")
            {
                double tmp=0;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _max_subjet_mass = tmp;
            }
            else if (key == "mode")
            {
                int tmp=4;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                set_mode(tmp);
            }
            else if (key == "mtmin" || key == "mintopmass")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _mtmin = tmp;
            }
            else if (key == "mtmax" || key == "maxtopmass")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _mtmax = tmp;
            }
            else if (key == "rmin")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _rmin = tmp;
            }
            else if (key == "rmax")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _rmax = tmp;
            }
            else if (key == "fw")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _fw = tmp;
                _rmin = (1.-_fw)*_mwmass/_mtmass;
                _rmax = (1.+_fw)*_mwmass/_mtmass;
            }
            else if (key == "m23cut" || key == "cut::m23")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _m23cut = tmp;
            }
            else if (key == "m13cutmin" || key == "cut::m13min")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _m13cutmin = tmp;
            }
            else if (key == "m13cutmax" || key == "cut::m13max")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _m13cutmax = tmp;
            }
            else if (key == "mintoppt" || key == "cut::ptmintop")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _minpt_tag = tmp;
            }
            else if (key == "nfiltering" || key == "filteringn" || key == "nfilt")
            {
                unsigned tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _filtering_n = tmp;
            }
            else if (key == "rfiltering" || key == "filteringr" || key == "rfilt")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _filtering_R = tmp;
            }
            else if (key == "filtering_minpt_subjet" || key == "filteringminptsubjet" || key == "minptfilt")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _filtering_minpT_subjet = tmp;
            }
            else if (key == "filtering_jetalgorithm" || key == "filteringjetalgorithm" || key == "jetalgfilt")
            {
                if (it->second == "antikT" || it->second == "antikt")
                {
                    _filtering_jetalg = fastjet::antikt_algorithm;
                }
                else if (it->second == "kT" || it->second == "kt")
                {
                    _filtering_jetalg = fastjet::kt_algorithm;
                }
                else if (it->second == "cambridge" || it->second == "cambridge_achen")
                {
                    _filtering_jetalg = fastjet::cambridge_algorithm;
                }
                else if (it->second == "genkt" || it->second == "genkT")
                {
                    _filtering_jetalg = fastjet::genkt_algorithm;
                }
            }
            else if (key == "reclustering_jetalgorithm" || key == "reclusteringjetalgorithm" || key == "jetalgreclust")
            {
                if (it->second == "antikT" || it->second == "antikt")
                {
                    _reclustering_jetalg = fastjet::antikt_algorithm;
                }
                else if (it->second == "kT" || it->second == "kt")
                {
                    _reclustering_jetalg = fastjet::kt_algorithm;
                }
                else if (it->second == "cambridge" || it->second == "cambridge_achen")
                {
                    _reclustering_jetalg = fastjet::cambridge_algorithm;
                }
                else if (it->second == "genkt" || it->second == "genkT")
                {
                    _reclustering_jetalg = fastjet::genkt_algorithm;
                }
            }
            else if (key == "cut::z" || key == "zcut")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _zcut = tmp;
            }
            else if (key == "cut::r_factor" || key == "rcut_factor")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _rcut_factor = tmp;
            }
            else if (key == "max_fatjet_r" || key == "maxfatjetr")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _max_fatjet_R = tmp;
            }
            else if (key == "min_fatjet_r" || key == "minfatjetr")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _min_fatjet_R = tmp;
            }
            else if (key == "step_r" || key == "stepr")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _step_R = tmp;
            }
            else if (key == "optimalr_treshold" || key == "optimalrtreshold")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _optimalR_threshold = tmp;
            }
            else if (key == "r_filt_optimalr_calc")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _R_filt_optimalR_calc = tmp;
            }
            else if (key == "n_filt_optimaln_calc")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _N_filt_optimalR_calc = tmp;
            }
            else if (key == "r_min_exp_function")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _r_min_exp_function = tmp;
            }
            else if (key == "optimalr_mmin")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _optimalR_mmin = tmp;
            }
            else if (key == "optimalr_mmax")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _optimalR_mmax = tmp;
            }
            else if (key == "optimalr_fw")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _optimalR_fw = tmp;
            }
            else if (key == "r_opt_diff")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _R_opt_diff = tmp;
            }
            else if (key == "r_filt_optimalr_pass")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _R_filt_optimalR_pass = tmp;
            }
            else if (key == "n_filt_optimalr_pass")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _N_filt_optimalR_pass = tmp;
            }
            else if (key == "r_filt_optimalr_fail")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _R_filt_optimalR_fail = tmp;
            }
            else if (key == "n_filt_optimalr_fail")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _N_filt_optimalR_fail = tmp;
            }
            else if (key == "q_zcut" || key=="cut::q_z")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _q_zcut = tmp;
            }
            else if (key == "q_dcut_fctr" || key=="cut::q_d_fctr")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _q_dcut_fctr = tmp;
            }
            else if (key == "q_exp_min")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _q_exp_min = tmp;
            }
            else if (key == "q_exp_max")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _q_exp_max = tmp;
            }
            else if (key == "q_rigidity")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _q_rigidity = tmp;
            }
            else if (key == "q_truncation_fctr")
            {
                double tmp=0.;
                std::stringstream str;
                str << it->second;
                str >> tmp;
                _q_truncation_fctr = tmp;
            }
        }
}


}
