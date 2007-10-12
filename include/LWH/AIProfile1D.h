// -*- C++ -*-
#ifndef LWH_AIProfile1D_H
#define LWH_AIProfile1D_H 
//

#ifndef LWH_USING_AIDA

/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class IAnnotation;


class IBaseHistogram {

public:

  virtual ~IBaseHistogram() {}

  virtual std::string title() const = 0;
  virtual bool setTitle(const std::string & title) = 0;
  virtual int dimension() const = 0;
  virtual bool reset() = 0;
  virtual int entries() const = 0;

};

class IProfile : virtual public IBaseHistogram {

public:

    virtual ~IProfile() {}

    virtual int allEntries() const = 0;
    virtual int extraEntries() const = 0;
    virtual double sumBinHeights() const = 0;
    virtual double sumAllBinHeights() const = 0;
    virtual double sumExtraBinHeights() const = 0;
    virtual double minBinHeight() const = 0;
    virtual double maxBinHeight() const = 0;

};

class IAxis;

class IProfile1D : virtual public IProfile {

public: 

    virtual ~IProfile1D() {}

    virtual bool fill(double x, double y, double weight = 1.) = 0;
    virtual double binMean(int index) const = 0;
    virtual int binEntries(int index) const = 0;
    virtual double binHeight(int index) const = 0;
    virtual double binError(int index) const = 0;
    virtual double binRms(int index) const = 0;
    virtual double mean() const = 0;
    virtual double rms() const = 0;
    virtual const IAxis & axis() const = 0;
    virtual int coordToIndex(double coord) const = 0;
    virtual bool add(const IProfile1D & profile) = 0;

}; // class

} // namespace AIDA

/** @endcond */
  
#else
#include "AIDA/IProfile1D.h"
#endif

#endif /* LWH_AIProfile1D_H */

