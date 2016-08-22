// -*- C++ -*-
#ifndef RIVET_Histo1D_HH
#define RIVET_Histo1D_HH


#include "YODA/AnalysisObject.h"
#include "YODA/Histo1D.h"

namespace YODA {
    typedef std::shared_ptr<YODA::Histo1D> Histo1DPtr;
}

#include <HepMC/WeightContainer.h>

namespace Rivet {

    typedef shared_ptr<YODA::AnalysisObject> AnalysisObjectPtr;

    class Histo1D : public YODA::AnalysisObject {
        public:

            Histo1D(std::vector<std::string> _weightNames, const YODA::Scatter2D& refscatter, const string& path); 
            Histo1D(std::vector<std::string> _weightNames, const vector<double>& binedges, const string& path);
            Histo1D(std::vector<std::string> _weightNames, size_t nbins, double lower, double upper, const string& path);

            void fill(double x, double weight) {
                _active->fill(x,weight);
            }

            const YODA::Histo1D & activeobj() const { assert(_active); return *_active; }
            YODA::Histo1DPtr active() const { return _active; }

            double integral() { assert(false); return 0.0; }

            const vector<YODA::Histo1DPtr> & persistent() const { return _persistent; }

            size_t dim() const { return 1; }

        private:
            // @todo
            // is this the right thing to do?
            friend class AnalysisHandler;

            // ?? automatically ??
            //	void new_Active(YODA::Histo1D * newA) {
            //		_active = newA;
            //	}

            void newEvent() {
                YODA::Histo1DPtr tmp = make_shared<YODA::Histo1D>(_persistent[0]->clone());
                tmp->reset();
                _evgroup.push_back( tmp );
                _active = _evgroup.back();
            }


            void setActiveWeight(unsigned int iWeight) {
                _active = _persistent.at(iWeight);
            }


            void unsetActiveWeight() {
                _active.reset();
            }

            void pushToPersistent(const vector<vector<double> >& weight);

            void reset() { 
                foreach (AnalysisObjectPtr ao, _persistent) {
                    ao->reset();
                }
                _evgroup.clear();
                _active.reset(); 
            }


            AnalysisObject* newclone() const {
                return new Histo1D(*this);
            }

        private:

            vector<YODA::Histo1DPtr> _persistent; // M of these, one for each weight

            vector<YODA::Histo1DPtr> _evgroup; // N of these, one for each event in ev.group

            YODA::Histo1DPtr _active;

    };

    inline YODA::Scatter2D toIntegralHisto(const Histo1D& h, bool includeunderflow=true) {
        return toIntegralHisto(h.activeobj(),includeunderflow);
    }


    inline YODA::Scatter2D divide(const Histo1D& numer, const Histo1D& denom) {
        return divide(numer.activeobj(),denom.activeobj());
    }

    inline YODA::Scatter2D operator / (const Histo1D& numer, const Histo1D& denom) {
        return divide(numer,denom);
    }

}

#endif
