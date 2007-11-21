// -*- C++ -*-
#ifndef RIVET_SVertex_HH
#define RIVET_SVertex_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Event.hh"


namespace Rivet {


  /**
     @brief Determine secondary vertices.
     
     Makes use of PVertex projection.

     @todo Replace function with a functor to improve equality comparisons.

     Complex cuts on tracks and vertices to validate them have to be provided
     by an external function
     bool f(SVertex&, ParticleVector&, const HepMC::GenVertex&, FourMomentum);
     which can be embedded in the analysis code. An example can be found 
     in the HepEx0605099 analysis. A pointer to this function has to be given 
     to the constructor of the SVertex projection. Its arguments are as follows:

     in: reference to instance of SVertex projection, ParticleVector of
         vertex to be analyzed, primary (Gen)Vertex
     out: FourMomentum = visible Momentum of vertex (selected tracks), 
     return bool: cuts passed? 1 : 0 

     In this way the SVertex projection can be kept as universal/flexible
     as possible.

     The constructor expects also a list of (pre-selected) jets.
     Associated tracks and vertices to a jet are checked for displacement.
     A list of tagged jets can be obtained via the getTaggedJets() function
  */

  class SVertex: public PVertex {

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor. Must specify a PVertex 
    /// projection object which is assumed to live through the run.
    SVertex(PVertex& pvtx, ChargedFinalState& chfs, 
            vector<FourMomentum>& jetaxes, double deltaR,
            bool (*applyVtxTrackCuts) (SVertex&, ParticleVector&, const HepMC::GenVertex&, FourMomentum),
            double detEta, 
            double IPres, double DLS, 
            double DLSres=0.) 
      : _pvtx(pvtx), _chfs(chfs), _jetaxes(jetaxes), _deltaR(deltaR),
        _applyVtxTrackCuts(applyVtxTrackCuts),
        _detEta(detEta), _IPres(IPres), _DLS(DLS), 
        _DLSres(DLSres)
    { 
      addProjection(_pvtx);
      addProjection(_chfs);
      if (_DLSres==0.) _DLSres=_IPres;
    }

    /// The destructor.
    virtual ~SVertex() { }
    //@}

  public:
    /// Return the name of the projection
    string getName() const {
      return "SVertex";
    }

    /// Return vector of tagged jets (FourMomentum's)
    const vector<FourMomentum>& getTaggedJets() const {
      return _taggedjets;
    }

    /// Return Distance of Closest Approach from track to given (primary) vertex
    double get2dDCA(const HepMC::GenParticle& track, const HepMC::GenVertex& gvtx);

    /// Return Impact Parameter Significance of givern track w.r.t. (primary) vertex
    double get2dDCAsig(const HepMC::GenParticle& track, const HepMC::GenVertex& gvtx);

    /// Return Distance of Closest Approach from track to given (primary) vertex
    double get3dDCA(const HepMC::GenParticle& track, const HepMC::GenVertex& gvtx);

    /// Return Impact Parameter Significance of givern track w.r.t. (primary) vertex
    double get3dDCAsig(const HepMC::GenParticle& track, const HepMC::GenVertex& gvtx);

    /// Return Decay Length Significance between two vertices in transverse plane
    double get2dDLS(const HepMC::GenVertex& vtx1, const HepMC::GenVertex& vtx2, 
		    FourMomentum& jetaxis);

    /// Return 3 dimensional Decay Length Significance between vertices 
    double get3dDLS(const HepMC::GenVertex& vtx1, const HepMC::GenVertex& vtx2, 
		    FourMomentum& jetaxis);


  protected:

    /// Apply the projection to the event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  private:

    /// The Primary Vertex projection used by this projection
    PVertex _pvtx;

    /// The ChargedFinalState projection used by this projection
    ChargedFinalState _chfs;

    /// The jet axes of the jet algorithm projection
    vector<FourMomentum> _jetaxes;

    /// Max distance between vis. momentum of vertex and jet to be probed
    double _deltaR;

    /// Analysis dependent cuts to be specified in analysis function
    bool (*_applyVtxTrackCuts) (SVertex&, ParticleVector&, const HepMC::GenVertex&, FourMomentum);

    /// Geometrical acceptance of tracker
    double _detEta;

    /// Impact parameter resolution, (including beam size)
    double _IPres;

    /// Decay length significance (cut value)
    double _DLS;

    /// Decay length significance uncertainty
    double _DLSres;

    /// Visible Momentum 4-vector of given vertex
    //FourMomentum _vtxVisMom;

    /// Jets which have been tagged
    vector<FourMomentum> _taggedjets;
  };

}

#endif
