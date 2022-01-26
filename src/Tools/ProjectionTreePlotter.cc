// -*- C++ -*-
//ProjectionTreePlotter
//Utility to plot the "tree" of projections of a given analysis.

#include "Rivet/Tools/ProjectionTreePlotter.hh"

#include <fstream>

//Small class used to store nodes of the projection tree.
//OOP is possibly a bit of overkill but it works and nothing else needs to know about this code, so its fine I hope.
class ProjectionTreeNode{
private:
    //A node has EITHER an AnaHandle or ProjectionPtr.
    Rivet::ConstProjectionPtr _theptr;
    Rivet::AnaHandle _anahandle;
    //Has this node had all its children examined
    mutable bool _analysed;
    //At what index is this node stored in the list of nodes.
    const size_t _location;
public:
    //Standard constructor
    ProjectionTreeNode(size_t loc, Rivet::ConstProjectionPtr ProjPointer, Rivet::AnaHandle Anapointer = nullptr) : _location(loc){
        _theptr = ProjPointer;
        _anahandle = Anapointer;
        _analysed=false;
    }

    //getter function
    bool is_analysed() const {return _analysed;}

    //Get children from the Projection Handler
    std::set<Rivet::ConstProjectionPtr> getChildren() const {
        if (_theptr != nullptr) {return _theptr->getImmediateChildProjections();}
        else {return _anahandle->getImmediateChildProjections();}
    }

    ///Get projection/analysis's name
    std::string get_name() const {
        if (_theptr == nullptr){return _anahandle->name();}
        else {return _theptr->name();}
    }

    //Function used in building the tree.
    //Get children: if they've not already been analysed, push them to the end of the passed-by-reference vector of all projections.
    void add_children_to_list(std::vector<Rivet::ConstProjectionPtr>& MainProjVector, std::vector<std::array<size_t,2>>& MainEdgeVector, 
                                std::vector<std::string>& nameVector, std::vector<ProjectionTreeNode>& MainNodeVector) const {
        auto children = getChildren();
        if (children.size() == 0) {_analysed = true; return;}
        int n_analysed = 0;
        for (Rivet::ConstProjectionPtr child : children){
            auto it = std::find(MainProjVector.begin(), MainProjVector.end(), child);
            // If element was found
            if (it != MainProjVector.end())
            {
                size_t index = it - MainProjVector.begin();
                MainEdgeVector.push_back({_location, index});
                if (MainNodeVector[index].is_analysed()) {++n_analysed;}
            }
            else {
                // If the element is not present in the vector
                size_t index = MainProjVector.size();
                MainProjVector.push_back(child);
                MainEdgeVector.push_back({this->_location, index});
                ProjectionTreeNode newnode(index, child);
                MainNodeVector.push_back(newnode);
                nameVector.push_back(newnode.get_name());
            }
        }
        if (n_analysed == children.size()) {
            _analysed=true;
            return;
        }
    }

    //Get Logger object
    Rivet::Log& getLog() const {
        return Rivet::Log::getLog("Rivet.ProjectionTreeNode");
    }
};



namespace Rivet{

    ProjectionTreeGenerator::ProjectionTreeGenerator(const std::string& name) : _path(name), _analyses({}), _edgeVector({}), 
                                                                                _nameVector({}), _projVector({}), _title("") {
    }

    ProjectionTreeGenerator::ProjectionTreeGenerator() : ProjectionTreeGenerator("") {
    }

    void ProjectionTreeGenerator::setPath(const std::string&path){
        _path = path;
    }

    void ProjectionTreeGenerator::setTitle(const std::string&title){
        _title = title;
    }

    void ProjectionTreeGenerator::addAnalysis(const std::string& analysis){
        _analyses.push_back(analysis);
        _treeGenerated = false;
    }

    void ProjectionTreeGenerator::addAnalyses(const std::vector<std::string>& analyses){
        for (const std::string& ana : analyses){
            _analyses.push_back(ana);
        }
        _treeGenerated=false;
    }

    int ProjectionTreeGenerator::generateProjTree() {
        //The analysishandler for this ProjectionTree.
        Rivet::AnalysisHandler ah;
        ah.setCheckBeams(false);
        for (const std::string& analysis : _analyses){
            ah.addAnalysis(analysis);
        }

        // Create a dummy event to initialise with
        GenEvent e;
        ah.analyze(e);

        //Use the projectionTreeNode
        std::vector<ProjectionTreeNode> nodeVector;

        for (size_t i=0; i < ah.analyses().size(); ++i){
            _projVector.push_back(nullptr);
            ProjectionTreeNode ananode(i,nullptr,ah.analyses()[i]);
            nodeVector.push_back(ananode);
            _nameVector.push_back(ananode.get_name());
        }

        size_t num_analyses = nodeVector.size();
        MSG_INFO("Constructing Projection Tree for " << num_analyses << " analys" << ((num_analyses==1)?"is":"es"));

        //TODO: this would be more elegant with iterators
        size_t tracker = 0;
        do {
            nodeVector[tracker].add_children_to_list(_projVector, _edgeVector, _nameVector, nodeVector);
            tracker +=1;
        } while (tracker < nodeVector.size());

        MSG_INFO("Number of unique projections " << _projVector.size() - num_analyses);
        _treeGenerated = true;
        return 0;
    }


    void ProjectionTreeGenerator::saveProjTree() const {
        if (!_treeGenerated){
            MSG_WARNING("Trying to save a projection tree that has not yet been generated. Please check your code calls generate first!");
            return;
        }
        std::ofstream outfile;
        std::string digraphname = _title;
        //If no title set get from path.
        if (digraphname == ""){
            //Remove anything before the last '\'. There's probably a nice regex way but this is still only two lines.
            auto beginning = std::find(_path.rbegin(), _path.rend(), '/');
            digraphname = std::string(beginning.base(), _path.end());
            //Strip the .gv if it's there
            if (std::string(digraphname.end()-3,digraphname.end()) == ".gv"){
                digraphname = std::string(digraphname.begin(), digraphname.end()-3);
            }
        }

        outfile.open(_path);
        outfile << "digraph \"" << digraphname <<  "\"{\n";
        for (size_t i = 0; i < _nameVector.size(); ++i){
            //Do analysis boxes in a different colour:
            if (i < _analyses.size()){
                outfile << i << "[fillcolor=\"#F09C9C\", style=\"rounded,filled\", shape=box,label=<"<<_nameVector[i]<<">];\n";
            } else {
                outfile << i << "[fillcolor=\"#F0F0D0\", style=\"rounded,filled\", shape=box,label=<"<<_nameVector[i]<<">];\n";
            }
        }
        for (const auto & edge : _edgeVector){
            outfile << edge[0] << "->" << edge[1] << " ;\n";
        }
        outfile << "}";
        outfile.close();

        MSG_INFO("Saved Projection Tree to " << _path << " (use graphviz's dot or similar to produce an image file - e.g. \"dot -Tsvg "<<_path<<" > "<<digraphname<<".svg \")");
    }
    
    Log& ProjectionTreeGenerator::getLog() const {
        return Rivet::Log::getLog("Rivet.ProjectionTreeGenerator");
    }
}