// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/ProjectionHandler.hh"

#ifndef RIVET_RivetProjectionTree_HH
#define RIVET_RivetProjectionTree_HH

namespace Rivet{

  /// @brief Class that deals with generating projection trees (for debugging etc.)
  class ProjectionTreeGenerator{
    private:
      ///@name Private member variables that define what tree we build:
      ///{
      ///Which analyses to generate projection tree for
      std::vector<std::string> _analyses;

      ///path to save tree to.
      std::string _path;

      //Title of digraph (by default the last bit of the name minus the '.gv')
      std::string _title;

      ///flag to keep track of whether or not a tree has been generated
      bool _treeGenerated;
      ///}

      ///@name The vectors that store the projection tree itself.
      ///{
      ///Stores all the projections in the tree
      std::vector<Rivet::ConstProjectionPtr> _projVector;

      ///Stores all the edges: an edge is directional and stored as array<size_t, 2>{start, finish};
      std::vector<std::array<size_t, 2>> _edgeVector;

      ///Stores the name labels of the projections in _projVector
      std::vector<std::string> _nameVector; //Should be 1<->1 with _projVector;
      ///}
      

    public:
      ///Standard constructor with name of gv pre-supplied
      ProjectionTreeGenerator(const std::string& name);

      ///Standard constructor
      ProjectionTreeGenerator();

      //Set the path
      void setPath(const std::string& path);

      //Set the title (defaults to last bit of path without '.gv')
      void setTitle(const std::string& title);

      ///Add analysis
      void addAnalysis(const std::string& analysis);

      ///Add analyses
      void addAnalyses(const std::vector<std::string>& analyses);
    
      ///Generate the projection tree
      int generateProjTree();

      ///Save the projection tree to it's _path file
      void saveProjTree() const;

      /// Get a logger object.
      Log& getLog() const;
  };
}

#endif