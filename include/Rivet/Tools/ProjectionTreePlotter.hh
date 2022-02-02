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
      ///@name Private member variables that dictate output of the tree
      ///{
      ///path to save tree to.
      std::string _path;

      //Title of digraph (by default the last bit of the name minus the '.gv')
      std::string _title;

      ///flag to keep track of whether or not a tree has been generated
      bool _treeGenerated;

      //Tracks how many analyses are involved for the projection tree
      size_t _nAnalyses;
      ///}

      ///@name The vectors that store the projection tree itself.
      ///{
      ///Stores all the projections in the tree
      std::vector<Rivet::ConstProjectionPtr> _projVector;

      ///Stores all the edges: an edge is directional and stored as pair<size_t, size_t>{start, finish};
      ///(switched from array<size_t,2> to pair as cython doesn't wrap std::array)
      std::vector<std::pair<size_t, size_t>> _edgeVector;

      ///Stores the name labels of the projections in _projVector
      std::vector<std::string> _nameVector; //Should be 1<->1 with _projVector;
      ///}
      

    public:
      ///Standard constructor with name of gv pre-supplied
      ProjectionTreeGenerator(const std::string& path);

      ///Standard constructor
      ProjectionTreeGenerator();

      ///Set the path
      void setPath(const std::string& path);

      ///Set the title (defaults to last bit of path without '.gv')
      void setTitle(const std::string& title);
    
      ///Generate the projection tree for the supplied analyses (constructs dummy analysisHandler for you)
      int generateProjTree(const std::vector<std::string>& analyses);

      ///Get a projection tree from the supplied analysishandler
      int generateProjTree(const AnalysisHandler& ah);

      ///Save the projection tree to the path specifed by _path
      void writeProjTree() const;

      ///Get the vector of projection names
      ///TODO: I'd have preferred a pass-by-reference solution but cython wasn't co-operating.
      inline const std::vector<string>& getProjNames() const{
        return _nameVector;
      }

      ///Get the vector of edges - format pair<size_t,size_t>(start-index,end-index)
      inline const std::vector<std::pair<size_t, size_t>>& getEdges() const {
        return _edgeVector;
      }

      /// Get a logger object.
      Log& getLog() const;
  };
}

#endif