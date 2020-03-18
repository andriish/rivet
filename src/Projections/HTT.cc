// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/HTT.hh"

namespace Rivet {

using namespace std;


HTT::HTT(const JetAlg& jetalg)
{
    setName("HEPTopTagger");
    declare(jetalg, "Jets");
    set_mode(4);
}

CmpState HTT::compare(const Projection& p) const 
{

    return CmpState::EQ;
}

void HTT::Reset()
{
    _topjets.clear();
}

void HTT::calc(const Jets& jets) {

    for (unsigned i=0; i<jets.size();i++) {
      // Apply jet cuts
      HEPTopTagger::HEPTopTagger tagger(jets[i]);
             // Unclustering, Filtering & Subjet Settings
      tagger.set_max_subjet_mass(_max_subjet_mass);
      tagger.set_mass_drop_threshold(_mass_drop_treshold);
      tagger.set_filtering_R(_filtering_R);
      tagger.set_filtering_n(_filtering_n);
      tagger.set_filtering_minpt_subjet(_filtering_minpT_subjet);

      // How to select among candidates
      tagger.set_mode(_mode);

      // Requirements to accept a candidate
      tagger.set_top_minpt(_minpt_tag);
      tagger.set_top_mass_range(_mtmin, _mtmax);
      tagger.set_fw(_fw);

      // Run the tagger
      tagger.get_setting();
      tagger.get_info();
      tagger.run();
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

void HTT::project(const Event& e) {
    const Jets jets = applyProjection<JetAlg>(e, "Jets").jets();
    calc(jets);
}



//void HTT::set_mtmass(double v) {_mtmass=v;}
//void HTT::set_mwmass(double v) {_mwmass=v;}
//void HTT::set_mtmin(double v) {_mtmin=v;}
//void HTT::set_mtmax(double v) {_mtmax=v;}
//void HTT::set_rmin(double v) {_rmin=v;}
//void HTT::set_rmax(double v) {_rmax=v;}
//void HTT::set_fw(double fw) {_fw=fw;}
//
//void HTT::set_m23cut(double v) {_m23cut=v;}
//void HTT::set_m13cutmin(double v) {_m13cutmin=v;}
//void HTT::set_m13cutmax(double v) {_m13cutmin=v;}
//void HTT::set_minpt_tag(double v) {_minpt_tag=v;}
//
//void HTT::set_filtering_n(unsigned i) {_filtering_n=i;}
//void HTT::set_filtering_R(double v) {_filtering_R=v;}
//void HTT::set_filtering_minpT_subjet(double v) {_filtering_minpT_subjet=v;}
//
//void HTT::set_filtering_jetalg(fastjet::JetAlgorithm j) {_filtering_jetalg=j;}
//void HTT::set_reclustering_jetalg(fastjet::JetAlgorithm j) {_reclustering_jetalg=j;}



}