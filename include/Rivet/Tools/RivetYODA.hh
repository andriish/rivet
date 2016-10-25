#ifndef RIVET_RIVETYODA_HH
#define RIVET_RIVETYODA_HH

/// @author Andy Buckley
/// @date   2009-01-30
/// @author David Grellscheid
/// @date   2011-07-18
/// @author David Grellscheid
/// @date   2016-09-27

#include "Rivet/Config/RivetCommon.hh"

#include "YODA/AnalysisObject.h"
#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/Scatter1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Scatter3D.h"

#include <map>

namespace YODA {
    typedef std::shared_ptr<YODA::AnalysisObject> AnalysisObjectPtr;
    typedef std::shared_ptr<YODA::Histo1D> Histo1DPtr;
    typedef std::shared_ptr<YODA::Histo2D> Histo2DPtr;
    typedef std::shared_ptr<YODA::Profile1D> Profile1DPtr;
    typedef std::shared_ptr<YODA::Profile2D> Profile2DPtr;
    typedef std::shared_ptr<YODA::Scatter1D> Scatter1DPtr;
    typedef std::shared_ptr<YODA::Scatter2D> Scatter2DPtr;
    typedef std::shared_ptr<YODA::Scatter3D> Scatter3DPtr;
    typedef std::shared_ptr<YODA::Counter> CounterPtr;
}



namespace Rivet {


    class AnalysisObjectPtr {
        public:
            virtual ~AnalysisObjectPtr() {}

            virtual YODA::AnalysisObject* operator->() = 0;
            virtual YODA::AnalysisObject* operator->() const = 0;
            virtual const YODA::AnalysisObject & operator*() const = 0;

            bool operator ==(const AnalysisObjectPtr& p) { return (this == &p); }

        protected:
            /// @todo do we need this?
            // virtual void reset() = 0;
    };

    /// @todo
    /// implement scatter1dptr and scatter2dptr here
    class Scatter1DPtr : public AnalysisObjectPtr {
        public:
            Scatter1DPtr() :
                _scatter(YODA::Scatter1DPtr()) { }

            Scatter1DPtr(const YODA::Scatter1DPtr& p) :
                _scatter(p) { }

            bool operator!() const { return !_scatter; }
            operator bool() const { return bool(_scatter); }

            YODA::Scatter1D* operator->() { return _scatter.get(); }

            YODA::Scatter1D* operator->() const { return _scatter.get(); }

            YODA::Scatter1D & operator*() { return *_scatter; }

            const YODA::Scatter1D & operator*() const { return *_scatter; }

        protected:
            YODA::Scatter1DPtr _scatter;
    };

    class Scatter2DPtr : public AnalysisObjectPtr {
        public:
            Scatter2DPtr(const YODA::Scatter2DPtr& p) :
                _scatter(p) { }

            Scatter2DPtr() :
                _scatter(YODA::Scatter2DPtr()) { }

            bool operator!() { return !_scatter; }
            operator bool() { return bool(_scatter); }

            YODA::Scatter2D* operator->() { return _scatter.get(); }

            YODA::Scatter2D* operator->() const { return _scatter.get(); }

            YODA::Scatter2D & operator*() { return *_scatter; }

            const YODA::Scatter2D & operator*() const { return *_scatter; }

        protected:
            YODA::Scatter2DPtr _scatter;
    };

    class Scatter3DPtr : public AnalysisObjectPtr {
        public:
            Scatter3DPtr(const YODA::Scatter3DPtr& p) :
                _scatter(p) { }

            Scatter3DPtr() :
                _scatter(YODA::Scatter3DPtr()) { }

            bool operator!() { return !_scatter; }
            operator bool() { return bool(_scatter); }

            YODA::Scatter3D* operator->() { return _scatter.get(); }

            YODA::Scatter3D* operator->() const { return _scatter.get(); }

            YODA::Scatter3D & operator*() { return *_scatter; }

            const YODA::Scatter3D & operator*() const { return *_scatter; }

        protected:
            YODA::Scatter3DPtr _scatter;
    };


    class MultiweightAOPtr : public AnalysisObjectPtr {

        public:
            virtual void newSubEvent() = 0;

            /// @todo 
            /// rename to setActive(Idx)?
            virtual void setActiveWeightIdx(unsigned int iWeight) = 0;

            virtual void pushToPersistent(const vector<vector<double> >& weight) = 0;

            virtual YODA::AnalysisObjectPtr activeYODAPtr() const = 0;
    };



template <class T>
    using Fill = pair<typename T::FillType, double>;

template <class T>
    using Fills = multiset<Fill<T>>;


class Histo1DTuple : public YODA::Histo1D {
public:
    Histo1DTuple(const YODA::Histo1D & h) : YODA::Histo1D(h) {}

    // todo: do we need to deal with users using fractions directly?
    void fill(double x, double weight=1.0, double fraction=1.0) {
        fills_.insert( {x, weight} );
    }

    void reset() {
        fills_.clear();
    }

    const Fills<YODA::Histo1D> & fills() const { return fills_; }

private:
    // x / weight pairs 
    Fills<YODA::Histo1D> fills_;
};









#define RIVETAOPTR_COMMON(YODATYPE)                                            \
    typedef shared_ptr<YODATYPE##Tuple> YODATYPE##TuplePtr;                    \
                                                                               \
    class YODATYPE##Ptr : public MultiweightAOPtr {                            \
        public:                                                                \
                                                                               \
        /* @todo                                                               \
         * some things are not really well-defined here                        \
         * for instance: fill() in the finalize() method and integral() in     \
         * the analyze() method.                                               \
         */                                                                    \
                                                                               \
        YODATYPE##Ptr(size_t len_of_weightvec, const YODA::YODATYPE& p) {      \
            for (size_t i = 0; i < len_of_weightvec; i++)                      \
                _persistent.push_back(make_shared<YODA::YODATYPE>(p));         \
                                                                               \
            return;                                                            \
        }                                                                      \
                                                                               \
        YODA::YODATYPE##Ptr active() const { return _active; }                 \
                                                                               \
        /* @todo this probably need to loop over all? */                       \
        bool operator!() const {return !active();}                             \
        operator bool() const {return bool(active());}                         \
                                                                               \
        YODA::YODATYPE* operator->() {                                         \
            return active().get();                                             \
        }                                                                      \
                                                                               \
        YODA::YODATYPE* operator->() const {                                   \
            return active().get();                                             \
        }                                                                      \
                                                                               \
        YODA::YODATYPE & operator*() {                                         \
            return *active();                                                  \
        }                                                                      \
                                                                               \
        const YODA::YODATYPE & operator*() const {                             \
            return *active();                                                  \
        }                                                                      \
                                                                               \
        /* @todo                                                               \
         * these need to be re-thought out.                                    \
                                                                               \
         void reset() { active()->reset(); }                                   \
         */                                                                    \
                                                                               \
        /* @todo                                                               \
         * these probably need to loop over all?                               \
         * do we even want to provide equality?                                \
         */                                                                    \
        /* @todo                                                               \
         * how about no.                                                       \
        friend bool operator==(YODATYPE##Ptr a, YODATYPE##Ptr b){              \
            if (a._persistent.size() != b._persistent.size())                  \
                return false;                                                  \
                                                                               \
            for (size_t i = 0; i < a._persistent.size(); i++) {                \
                if (a._persistent.at(i) != b._persistent.at(i)) {              \
                    return false;                                              \
                }                                                              \
            }                                                                  \
                                                                               \
            return true;                                                       \
        }                                                                      \
                                                                               \
        friend bool operator!=(YODATYPE##Ptr a, YODATYPE##Ptr b){              \
            return !(a == b);                                                  \
        }                                                                      \
                                                                               \
                                                                               \
        friend bool operator<(YODATYPE##Ptr a, YODATYPE##Ptr b){               \
            if (a._persistent.size() >= b._persistent.size())                  \
                return false;                                                  \
                                                                               \
            for (size_t i = 0; i < a._persistent.size(); i++) {                \
                if (*(a._persistent.at(i)) >= *(b._persistent.at(i))) {        \
                    return false;                                              \
                }                                                              \
            }                                                                  \
                                                                               \
            return true;                                                       \
        }                                                                      \
        */                                                                     \
                                                                               \
                                                                               \
        private:                                                               \
        void setActiveWeightIdx(unsigned int iWeight) {                        \
            _active = _persistent.at(iWeight);                                 \
            return;                                                            \
        }                                                                      \
                                                                               \
        /* this is for dev only---we shouldn't need this in real runs. */      \
        void unsetActiveWeight() {                                             \
            _active.reset();                                                   \
            return;                                                            \
        }                                                                      \
                                                                               \
        void newSubEvent() {                                                   \
            YODATYPE##TuplePtr tmp                                             \
                = make_shared<YODATYPE##Tuple>(_persistent[0]->clone());       \
            tmp->reset();                                                      \
            _evgroup.push_back( tmp );                                         \
            _active = _evgroup.back();                                         \
            return;                                                            \
        }                                                                      \
                                                                               \
        virtual YODA::AnalysisObjectPtr activeYODAPtr() const {                \
            return _active;                                                    \
        }                                                                      \
                                                                               \
        const vector<YODA::YODATYPE##Ptr> & persistent() const {               \
            return _persistent;                                                \
        }                                                                      \
                                                                               \
        /* to be implemented for each type */                                  \
        void pushToPersistent(const vector<vector<double> >& weight);          \
                                                                               \
        /* M of these, one for each weight */                                  \
        vector<YODA::YODATYPE##Ptr> _persistent;                               \
                                                                               \
        /* N of these, one for each event in evgroup */                        \
        vector<YODATYPE##TuplePtr> _evgroup;                                   \
                                                                               \
        YODA::YODATYPE##Ptr _active;                                           \
                                                                               \
        friend class AnalysisHandler;                                          \
    };

    // every object listed here needs a virtual fill method in YODA,
    // otherwise the Tuple fakery won't work.

    RIVETAOPTR_COMMON(Histo1D)
    // RIVETAOPTR_COMMON(Histo2D)
    // RIVETAOPTR_COMMON(Profile1D)
    // RIVETAOPTR_COMMON(Profile2D)
    // RIVETAOPTR_COMMON(Counter)

    using YODA::Counter;
    using YODA::Histo1D;
    using YODA::HistoBin1D;
    using YODA::Histo2D;
    using YODA::HistoBin2D;
    using YODA::Profile1D;
    using YODA::ProfileBin1D;
    using YODA::Profile2D;
    using YODA::ProfileBin2D;
    using YODA::Scatter1D;
    using YODA::Point1D;
    using YODA::Scatter2D;
    using YODA::Point2D;
    using YODA::Scatter3D;
    using YODA::Point3D;

    /// Function to get a map of all the refdata in a paper with the
    /// given @a papername.
    map<string, YODA::AnalysisObjectPtr> getRefData(const string& papername);

    /// @todo Also provide a Scatter3D getRefData() version?

    /// Get the file system path to the reference file for this paper.
    string getDatafilePath(const string& papername);

    /// Return the integral over the histogram bins
    /// @deprecated Prefer to directly use the histo's integral() method.
    DEPRECATED("Prefer to directly use the histo's integral() method.")
        inline double integral(Histo1DPtr histo) {
            return histo->integral();
        }


}

#endif
