#ifndef HepEntity_class
#define HepEntity_class

#include "inline_maths.h"


class HepEntity {

 public:

  HepEntity() {
    E=0.;
    px=0.;
    py=0.;
    pz=0.;
    /*
    pt=0.;
    p=0.;
    azi=0.;
    rap=0.;
    */
    return;
  }


  HepEntity(double E_in, double px_in, double py_in, double pz_in) : 
    E(E_in), px(px_in), py(py_in), pz(pz_in) {
    /*
    pt = sqrt(sqr(px)+sqr(py));
    p = sqrt(sqr(pt)+sqr(pz));
    azi = inline_maths::phi(px,py);
    rap = inline_maths::y(E,pz);
    */
    return;
  }


  //HepEntity(const HepEntity& in) : E(in.E), px(in.px), py(in.py), pz(in.pz),
  //pt(in.pt), p(in.p), azi(in.azi), rap(in.rap) {
  HepEntity(const HepEntity& in) : E(in.E), px(in.px), py(in.py), pz(in.pz) {
    return;
  }

  
  inline double y() const {
    //rap = inline_maths::y(E,pz);
    //return rap;
    return inline_maths::y(E,pz);
  }


  inline double phi() const {
    //azi = inline_maths::phi(px,py);
    //return azi;
    return inline_maths::phi(px,py);
  }

  inline double pT() const {
    //pt = sqrt(sqr(px)+sqr(py));
    //return pt;
    return sqrt(sqr(px)+sqr(py));
  }

  /*  
  inline double p() const {
    return sqrt(sqr(px)+sqr(py)+sqr(pz));
  }
  
  inline double m() const {
    double m2 = (sqr(E)-sqr(px)-sqr(py)-sqr(pz));
    return (m2>0) ? m2 : 0.;
  }
  */


  inline void p4vec(float* p) const {
    p[0] = px;
    p[1] = py;
    p[2] = pz;
    p[3] = E;
    return;
  }

  inline void Add(const HepEntity el) {
    E += el.E;
    px += el.px;
    py += el.py;
    pz += el.pz;
    /*
    pt = sqrt(sqr(px)+sqr(py));
    p = sqrt(sqr(pt)+sqr(pz));
    azi = inline_maths::phi(px,py);
    rap = inline_maths::y(E,p);
    */
    return;
  }

  inline void Fill(double E_in, double px_in, double py_in, double pz_in) {
    E = E_in;
    px = px_in;
    py = py_in;
    pz = pz_in;
    /*
    pt = sqrt(sqr(px)+sqr(py));
    p = sqrt(sqr(pt)+sqr(pz));
    azi = inline_maths::phi(px,py);
    rap = inline_maths::y(E,pz);
    */
    return;
  }


  double E;
  //double p;
  //double pt;
  double px;
  double py;
  double pz;
  //double rap;
  //double azi;

 private:



};
//end of class HepEntity;

#endif
