// -*- C++ -*-
#ifndef RIVET_EventMixingFinalState_HH
#define RIVET_EventMixingFinalState_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include <deque>
#include <algorithm>


namespace Rivet {

  // @brief Projects out an event mixed of several events, given
  // a mixing observable (eg. number of final state particles),
  // defining what should qualify as a mixable event.
  // Binning in the mixing observable is defined in the constructor,
  // as is the number of events one wants to mix with.
  // The method calculateMixingObs() must can be overloaded
  // in derived classes, to provide the definition of the mixing observables,
  // on the provided projection, eg. centrality or something more elaborate.
  //
  // !!!DISCLAIMER!!! The projection makes no attempt at correct handling
  // of event weights - ie. what one should use as event weight for several
  // mixed events. Proceed with caution if you do not use an MC with
  // unit weights.
  //
  // @author Christian Bierlich <christian.bierlich@thep.lu.se>

  typedef map<double, deque<Particles> > MixMap;	
  // EventMixingBase is the base class for event mixing projections.
  // Most methods are defined in this base class as they should.
  // In order to use it, a derived class should be implemented where:
  // - The constructor is reimplmented, giving the derived projection type
  // from which the mixing observable is calculated. The constructor must also
  // be declared public in the derived class.
  // - The calculateMixingObs is implemented.
  // To examples of such derived classes are given below,
  // 1) EventMixingFinalState, where the mixing observable are calculated
  // on a multiplicity of a charged final state, and:
  // 2) EventMixingCentrality, where the mixing observable is centrality.

  class EventMixingBase : public Projection {
  protected:
    // Constructor
    EventMixingBase(const Projection* mixObsProjPtr, const ParticleFinder& mix, size_t nMixIn,
	  double oMin, double oMax, double deltao ) : nMix(nMixIn){
    	setName("EventMixingBase");
	addProjection(*mixObsProjPtr,"OBS");
	addProjection(mix,"MIX");
	MSG_WARNING("EventMixing is a naive implementation, not currently " <<
	  "validated. Use with caution.");

	// Set up the map for mixing events.
	for(double o = oMin; o < oMax; o+=deltao )
	  mixEvents[o] = deque<Particles>();
	
    }

  public:
  
    // Return a vector of mixing events.
    vector<Particles> getMixingEvents() const {
	MixMap::const_iterator mixItr = mixEvents.lower_bound(mObs);
	if(mixItr == mixEvents.end() || mixItr->second.size() < nMix + 1)
	  return vector<Particles>();
	return vector<Particles>(mixItr->second.begin(), mixItr->second.end() - 1);
    }
 
  protected:

    // Calulate mixing observable.
    // Must be overloaded in derived classes.
    virtual void calculateMixingObs(const Projection* mProj) = 0;
    
    /// Perform the projection on the Event.
    void project(const Event& e){
	const Projection* mixObsProjPtr = &applyProjection<Projection>(e, "OBS");
	calculateMixingObs(mixObsProjPtr);
	MixMap::iterator mixItr = mixEvents.lower_bound(mObs);
	if(mixItr == mixEvents.end()){
	  // We are out of bounds.
	  MSG_DEBUG("Mixing observable out of bounds.");
	  return;
	}
	const Particles mix = applyProjection<ParticleFinder>(e, "MIX").particles();
	mixItr->second.push_back(mix);
	if(mixItr->second.size() > nMix + 1)
	  mixItr->second.pop_front();
    }

    /// Compare with other projections.
    int compare(const Projection& p) const {
	return mkNamedPCmp(p,"OBS");
    }
    
  // The mixing observable of the current event.
    double mObs;
  
  private:
    // The number of event to mix with.
    size_t nMix;
    // The event map.
    MixMap mixEvents;
  };

  // EventMixingFinalState has multiplicity in the mixing projection.
  // as the mixing observable.
  class EventMixingFinalState : public EventMixingBase {
   public:
    EventMixingFinalState(const ParticleFinder* mixObsProjPtr, const ParticleFinder& mix, size_t nMixIn,
      double oMin, double oMax, double deltao ) : 
      EventMixingBase(mixObsProjPtr, mix, nMixIn, oMin, oMax, deltao) {
    	setName("EventMixingFinalState");
      }

    DEFAULT_RIVET_PROJ_CLONE(EventMixingFinalState);
   
   protected:
    // Calulate mixing observable.
    virtual void calculateMixingObs(const Projection* mProj) {
      mObs = ((ParticleFinder*) mProj)->particles().size();
    }
  };

  // EventMixingCentrality has centrality as the mixing observable.
  class EventMixingCentrality : public EventMixingBase {
    public:
      EventMixingCentrality(const CentralityProjection* mixObsProjPtr, const ParticleFinder& mix, size_t nMixIn,
        double oMin, double oMax, double deltao ) : 
      EventMixingBase(mixObsProjPtr, mix, nMixIn, oMin, oMax, deltao) {
    	setName("EventMixingCentrality");
      }
    
      DEFAULT_RIVET_PROJ_CLONE(EventMixingCentrality);
    protected:
      virtual void calculateMixingObs(const Projection* mProj) {
        mObs = ((CentralityProjection*) mProj)->operator()();
      }
  }; 
}
#endif
