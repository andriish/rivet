// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/HTT.h"

using namespace Rivet;
using namespace std;


HTT::HTT(const JetAlg& jetalg)
{
    setName("HEPTopTagger");
    declare(jetalg, "Jets");
  _tagger.set_max_subjet_mass(30.);
  _tagger.set_mass_drop_threshold(0.8);
  _tagger.set_filtering_R(0.3);
  _tagger.set_filtering_n(5);
  _tagger.set_filtering_minpt_subjet(30.);

  // How to select among candidates
  _tagger.set_mode(HEPTopTagger::TWO_STEP_FILTER);

  // Requirements to accept a candidate
  _tagger.set_top_minpt(200);
  _tagger.set_top_mass_range(150., 200.);
  _tagger.set_fw(0.15);
  _tagger.set_debug(1);
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
    Reset();

    
//    for (const Jet& jet : jets) {
    for (unsigned i=0; i<jets.size();i++) {
      // Apply jet cuts
//      HEPTopTagger::HEPTopTagger tagger(jets[i]);
        _tagger(jets[i])
             // Unclustering, Filtering & Subjet Settings
//      tagger.set_max_subjet_mass(30.);
//      tagger.set_mass_drop_threshold(0.8);
//      tagger.set_filtering_R(0.3);
//      tagger.set_filtering_n(5);
//      tagger.set_filtering_minpt_subjet(30.);
//
//      // How to select among candidates
//      tagger.set_mode(HEPTopTagger::TWO_STEP_FILTER);
//
//      // Requirements to accept a candidate
//      tagger.set_top_minpt(200);
//      tagger.set_top_mass_range(150., 200.);
//      tagger.set_fw(0.15);

      // Run the tagger
      _tagger.run();
      MSG_INFO("Maybe top: " << _tagger.is_maybe_top());
            // Look at output if we have a tag:
      if (tagger.is_tagged()){
        MSG_INFO("Input fatjet: " << i << "  pT = " << jets[i].perp());
        MSG_INFO("Output: pT = " << _tagger.t().perp() << " Mass = " << _tagger.t().m() << " f_rec = " << _tagger.f_rec());
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