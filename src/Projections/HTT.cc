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

}