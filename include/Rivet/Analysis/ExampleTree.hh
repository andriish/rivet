// -*- C++ -*-
#ifndef RIVET_ExampleTree_H
#define RIVET_ExampleTree_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/KtJets.hh"
#include "Rivet/Projections/WZandh.hh"
#include "Rivet/RivetAIDA.fhh"

// Root stuff
#ifdef HAVE_ROOT
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#endif

namespace Rivet {

  /// This Analysis books and fills a ROOT tree with simulated data.
  /// Based initially on the ntuples used in Phys. Rev. D65; 096014 (2002)
  /// and JHEP05 (2007) 033.
  class ExampleTree : public Analysis {

  public:

    /// Default constructor
    inline ExampleTree()
      : p_fs(-4.0, 4.0, 0.0), p_chargedleptons(p_fs), p_ktjets(p_fs), p_wzandh()
    { 
      /// Particle IDs for neutrinos and antineutrinos and LSP
      p_vfs = new VetoedFinalState(p_fs);
      p_vfs->addVetoId(12,10.0,50.0);
      p_vfs->addVetoId(14);
      p_vfs->addVetoId(16);
      p_vfs->addVetoId(-12);
      p_vfs->addVetoId(-14);
      p_vfs->addVetoId(-16);
      p_vfs->addVetoId(1000022);
      p_totalvisiblemomentum = new TotalVisibleMomentum(*p_vfs);

      addProjection(p_fs);
      addProjection(p_chargedleptons);
      addProjection(p_ktjets);
      addProjection(*p_vfs);
      addProjection(*p_totalvisiblemomentum);
      addProjection(p_wzandh);
    }

    inline ~ExampleTree() {
      delete p_vfs;
      delete p_totalvisiblemomentum;
    }

  public:

    /// The name of this analysis is "Test"
    inline string getName() const {
      return "Test";
    }

  public:

    void init();
    
    void analyze(const Event & event);
    
    void finalize();

    /// Return the RivetInfo object of this analysis object.
    //    RivetInfo getInfo() const;

  private:

    /// The FinalState projector used by this analysis.
    FinalState p_fs;

    /// The Charged Lepton projector used by this analysis.
    ChargedLeptons p_chargedleptons;

    /// The jet projector
    KtJets p_ktjets;

    /// The vector boson projector
    WZandh p_wzandh;

    /// The VetoedFinalState projector used by this analysis.
    VetoedFinalState* p_vfs;

    /// The total visible momentum projector
    TotalVisibleMomentum* p_totalvisiblemomentum;



#ifdef HAVE_ROOT
    /// The tree
    TTree *rivetTree;

    /// The file for the Tree
    TFile *treeFile;

    /// The filename
    TString treeFileName;

#endif


    // The ntuple variables.
    int           nevt;            // event number

    int           nvb;             // number of W bosons
    float         vbvec[8][4];     // 4 momentum of W bosons.
    int           vbtype[8];       // type (i.e. decay mode) of W bosons.

    int           njet;            // number of jets
    float         vjet[50][4];     // four momentum of the jets

    int           nsub;            // number of jets for which the subjet analysis was performed.
    float         sjet3[200][4];   // four vector of jets for which we found subjets.
    float         ysubsj[200][4];     // y 1->2, 2->3, 3->4, 4->5 for the above jets.

    int           nlep;
    int           leptype[150][3];
    float         vlep[150][4];

    int           npart;           // Partons
    float         ppart[4000][4];
    int           pid[4000];
    int           mo[4000];

    float         esumr[4];        // Total visible momentum

  private:

    /// Hide the assignment operator
    ExampleTree & operator=(const ExampleTree& x);


  public:

    /// Minimum pt of jets which will go into the tree.
    int _jet_pt_cut;

    /// Minimum pt of jets which will have y evaluated and stored.
    int _subj_pt_cut;

    /// Minimum pt of charged leptons which will go into the tree.
    int _lepton_pt_cut;

    /// Store the partons or not?
    bool _store_partons;

  };

}

#endif
