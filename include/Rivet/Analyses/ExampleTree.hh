// -*- C++ -*-
#ifndef RIVET_ExampleTree_HH
#define RIVET_ExampleTree_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WZandh.hh"

// ROOT stuff
#ifdef HAVE_ROOT
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#endif


namespace Rivet {

  /// @brief Book and fill a ROOT tree with simulated data.
  ///
  /// This does some things, e.g. access parton level information, which
  /// are not recommended in rivet analyses, since the information is 
  /// unphysical and so cannot be compared to data, and also may be generator dependent.
  /// 
  class ExampleTree : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// Constructor
    ExampleTree();

    /// Factory method
    static Analysis* create() { return new ExampleTree(); }

    //@}


    /// @name Publication metadata
    //@{
    /// Return the name of this analysis
    string name() const {
      return "ExampleTree";
    }

    /// Get a description of the analysis.
    string spiresId() const {
      return "NONE";
    }

    /// Get a description of the analysis.
    string description() const {
      return "Demonstrate filling a ROOT tree from a kT jets analysis";
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "NONE";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "NONE";
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    #ifdef HAVE_ROOT
    /// The tree
    TTree* _rivetTree;
    
    /// The file for the Tree
    TFile* _treeFile;

    /// The filename
    TString _treeFileName;
    #endif
    

    /// @name The ntuple variables.
    //@{
    /// Event number
    int _nevt;            

    /// Number of W bosons
    int _nvb;             
    /// 4 momentum of W bosons.
    float _vbvec[8][4];
    /// Type (i.e. decay mode) of W bosons.
    int _vbtype[8]; 

    /// Number of jets
    int _njet; 
    /// Four momentum of the jets
    float _vjet[50][4]; 

    /// Number of jets for which the subjet analysis was performed.
    int _nsub; 
    /// Four vector of jets for which we found subjets.
    float _sjet3[200][4];
    /// y 1->2, 2->3, 3->4, 4->5 for the above jets.
    float _ysubsj[200][4];

    /// Number of leptons
    int _nlep;
    /// Lepton types
    int _leptype[150][3];
    float _vlep[150][4];

    /// Number of partons
    int _npart; 
    float _ppart[4000][4];
    int _pid[4000];
    int _mo[4000];

    /// Total visible momentum
    float _esumr[4];
    //@}


  private:

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
