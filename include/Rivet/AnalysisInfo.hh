// -*- C++ -*-
#ifndef RIVET_AnalysisInfo_HH
#define RIVET_AnalysisInfo_HH

#include "Rivet/Rivet.hh"
#include <ostream>

namespace Rivet {


  class AnalysisInfo {
    /// @todo Remove this when all metadata taken from YAML
    friend class Analysis;


  public:

    /// Static factory method: returns null pointer if no metadata found
    static AnalysisInfo* make(const std::string& name);

    /// @name Standard constructors and destructors.
    //@{

    /// The default constructor.
    AnalysisInfo() { clear(); }

    /// The destructor.
    ~AnalysisInfo() { }

    //@}


  public:

    /// @name Metadata
    /// Metadata is used for querying from the command line and also for
    /// building web pages and the analysis pages in the Rivet manual.
    //@{
    /// Get the name of the analysis. By default this is computed by
    /// combining the results of the experiment, year and Spires ID
    /// metadata methods and you should only override it if there's a
    /// good reason why those won't work.
    std::string name() const {
      if (!_name.empty()) return _name;
      if (!experiment().empty() && !year().empty() && !spiresId().empty()) {
        return experiment() + "_" + year() + "_S" + spiresId();
      }
      return "";
    }

    /// Get a description of the analysis.
    const std::string& spiresId() const { return _spiresId; }

    /// @brief Names & emails of paper/analysis authors.
    /// Names and email of authors in 'NAME <EMAIL>' format. The first
    /// name in the list should be the primary contact person.
    const std::vector<std::string>& authors() const { return _authors; }

    /// @brief Get a short description of the analysis.
    /// Short (one sentence) description used as an index entry.
    /// Use @a description() to provide full descriptive paragraphs
    /// of analysis details.
    const std::string& summary() const { return _summary; }

    /// @brief Get a full description of the analysis.
    /// Full textual description of this analysis, what it is useful for,
    /// what experimental techniques are applied, etc. Should be treated
    /// as a chunk of restructuredText (http://docutils.sourceforge.net/rst.html),
    /// with equations to be rendered as LaTeX with amsmath operators.
    const std::string& description() const { return _description; }

    /// @brief Information about the events needed as input for this analysis.
    /// Event types, energies, kinematic cuts, particles to be considered
    /// stable, etc. etc. Should be treated as a restructuredText bullet list
    /// (http://docutils.sourceforge.net/rst.html)
    const std::string& runInfo() const { return _runInfo; }

    /// Beam particle types
    const std::vector<std::pair<PdgId,PdgId> >& beams() const { return _beams; }

    /// Sets of valid beam energy pairs
    const std::vector<std::pair<double,double> >& energies() const { return _energies; }

    /// Experiment which performed and published this analysis.
    const std::string& experiment() const { return _experiment; }

    /// Collider on which the experiment ran.
    const std::string& collider() const { return _collider; }

    /// @brief When the original experimental analysis was published.
    /// When the refereed paper on which this is based was published,
    /// according to SPIRES.
    const std::string& year() const { return _year; }

    /// Journal, and preprint references.
    const std::vector<std::string>& references() const { return _references; }

    /// BibTeX citation key for this article.
    const std::string& bibKey() const { return _bibKey;}

    /// BibTeX citation entry for this article.
    const std::string& bibTeX() const {
      //return "@Article{" + bibKey() + ",\n" + _bibTeXBody + "\n}";
      return _bibTeX;
    }

    /// Whether this analysis is trusted (in any way!)
    const std::string& status() const { return _status; }

    /// Any work to be done on this analysis.
    const std::vector<std::string>& todos() const { return _todos; }
    //@}

    /// Return true if this analysis needs to know the process cross-section.
    bool needsCrossSection() const { return _needsCrossSection; }


  private:

    std::string _name;
    std::string _spiresId;
    std::vector<std::string> _authors;
    std::string _summary;
    std::string _description;
    std::string _runInfo;
    std::string _experiment;
    std::string _collider;
    std::vector<std::pair<PdgId, PdgId> > _beams;
    std::vector<std::pair<double, double> > _energies;
    std::string _year;
    std::vector<std::string> _references;
    std::string _bibKey;
    std::string _bibTeX;
    //std::string _bibTeXBody; //< Was thinking of avoiding duplication of BibKey...
    std::string _status;
    std::vector<std::string> _todos;
    bool _needsCrossSection;

    void clear() {
      _name = "";
      _spiresId = "";
      _authors.clear();
      _summary = "";
      _description = "";
      _runInfo = "";
      _experiment = "";
      _collider = "";
      _beams.clear();
      _energies.clear();
      _year = "";
      _references.clear();
      _bibKey = "";
      _bibTeX = "";
      //_bibTeXBody = "";
      _status = "";
      _todos.clear();
      _needsCrossSection = false;
    }

  };


  /// String representation
  std::string toString(const AnalysisInfo& ai);

  /// Stream an AnalysisInfo as a text description
  inline std::ostream& operator<<(std::ostream& os, const AnalysisInfo& ai) {
    os << toString(ai);
    return os;
  }


}

#endif
