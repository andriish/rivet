// -*- C++ -*-
#ifndef RIVET_DISRapidityGap_HH
#define RIVET_DISRapidityGap_HH

#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Get the incoming and outgoing hadron in a diffractive ep event
  class DISRapidityGap : public Projection {
  public:

    /// Type of DIS boost to apply
    enum Frame { HCM, LAB, XCM };

    DISRapidityGap() {
      setName("DISRapidityGap");
      declare(DISKinematics(), "DISKIN");
      declare(DISFinalState(DISFinalState::BoostFrame::HCM), "DISFS");
    }

    DEFAULT_RIVET_PROJ_CLONE(DISRapidityGap);

    /// @todo Document
    double t() const { return _t; }

    /// The absolute size of the largest rapidity gap
    double gap() const { return _gap; }

    /// The upper edge of the largest gap, relative to the DIS event orientation
    double gapUpper() const { return _gapUpp; }

    /// The lower edge of the largest gap, relative to the DIS event orientation
    double gapLower() const { return _gapLow; }

    /// @todo Document
    double epPzX(Frame f) const {
      if (f == LAB) return _ePpzX_LAB;
      else if (f == XCM) return _ePpzX_XCM;
      else return _ePpzX_HCM;
    }

    /// @todo Document
    const double emPzX(Frame f) const {
      if (f == LAB) return _eMpzX_LAB;
      else if (f == XCM) return _eMpzX_XCM;
      else return _eMpzX_HCM;
    }


    /// The particles defining system X
    ///
    /// @todo Document
    const Particles& systemX(Frame f) const {
      if (f == LAB) return _pX_LAB;
      else if (f == XCM) return _pX_XCM;
      else return _pX_HCM;
    }

    /// The particles defining system Y
    ///
    /// @todo Document
    const Particles& systemY(Frame f) const {
      if (f == LAB) return _pY_LAB;
      else if (f == XCM) return _pY_XCM;
      else return _pY_HCM;
    }

    /// Four-momentum of system X
    ///
    /// @todo Document
    const FourMomentum pX(Frame f) const {
      if (f == LAB) return _momX_LAB;
      else if (f == XCM) return _momX_XCM;
      else return _momX_HCM;
    }

    /// Four-momentum of system Y
    ///
    /// @todo Document
    const FourMomentum pY(Frame f) const {
      if (f == LAB) return _momY_LAB;
      else if (f == XCM) return _momY_XCM;
      else return _momY_HCM;
    }

    /// Mass^2 of system X
    ///
    /// @todo Document
    double m2X() const { return _M2X; }

    /// Mass^2 of system Y
    ///
    /// @todo Document
    double m2Y() const { return _M2Y; }


  protected:

    virtual CmpState compare(const Projection& p) const;

    virtual void project(const Event& e);

    void clear();

    void findGap(const Particles& particles, const DISKinematics& diskin);


  private:

    double _M2X, _M2Y, _t;
    double _gap, _gapUpp, _gapLow;
    double _ePpzX_LAB, _eMpzX_LAB;
    double _ePpzX_HCM, _eMpzX_HCM;
    double _ePpzX_XCM, _eMpzX_XCM;
    FourMomentum _momX_HCM, _momY_HCM;
    FourMomentum _momX_LAB, _momY_LAB;
    FourMomentum _momX_XCM, _momY_XCM;
    Particles _pX_HCM, _pY_HCM, _pX_LAB, _pY_LAB, _pX_XCM, _pY_XCM;

  };


}

#endif
