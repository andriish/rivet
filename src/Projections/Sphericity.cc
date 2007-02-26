// -*- C++ -*-

#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Cmp.hh"
#include "Rivet/RivetCLHEP.hh"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

//#include <Matrix.h>

using namespace Rivet;


int Sphericity::compare(const Projection & p) const {
  const Sphericity & other = dynamic_cast<const Sphericity &>(p);
  return pcmp(*fsproj, *other.fsproj);
}

RivetInfo Sphericity::getInfo() const {
  return Projection::getInfo() + fsproj->getInfo();
}

void Sphericity::project(const Event & e) {

  // Reset Parameters
  sphericity_ = 0;  
  planarity_ = 7;
  aplanarity_ = 0;

  const FinalState& fs = e.addProjection(*fsproj);
 
  CLHEP::HepMatrix mMom(3,3,0);
  //  map<int, const char*> cMap;
  //  cMap[1]="x()"; cMap[2]="y()"; cMap[3]="z()";
  double totalMomentumSq = 0;

  // iterate over all the final state particles
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
  
    // get the momentum vector for the final state particle
    CLHEP::HepLorentzVector lv = p->momentum;
 
    // if the original vector ia a HepVector instead
    // of lorentz then could do this:  
    // CLHEP::HepSymMatrix symMat =  hepv.T()*hepv;

    // total the momentum so that it can be used to 
    // normalise the momentum tensor
    totalMomentumSq += lv.x()*lv.x() + lv.y()*lv.y() + lv.z()*lv.z();

    // for now this crude method works at least 
    mMom[0][0]+=(lv.x()*lv.x()); 
    mMom[0][1]+=(lv.y()*lv.x()); 
    mMom[0][2]+=(lv.z()*lv.x()); 
    mMom[1][0]+=(lv.x()*lv.y());
    mMom[1][1]+=(lv.y()*lv.y());
    mMom[1][2]+=(lv.z()*lv.y());
    mMom[2][0]+=(lv.x()*lv.z());
    mMom[2][1]+=(lv.y()*lv.z()); 
    mMom[2][2]+=(lv.z()*lv.z());
  }

  mMom /= totalMomentumSq;

  // Check that the matrix is symmetric and if it is convert it,
  // diagonlize it, and rearrange the eigenvalues into the correct
  // order. Return the event shapes
  if   (mMom[0][1] == mMom[1][0]
    &&  mMom[0][2] == mMom[2][0]
    &&  mMom[1][2] == mMom[2][1]) {
    CLHEP::HepSymMatrix symMat;
    symMat.assign(mMom);
    CLHEP::diagonalize(&symMat);
    //std::cout << mMom << std::endl << std::endl;
    //std::cout << symMat << std::endl << std::endl;

    // Put the eigenvalues in the correct order
    vector<double> order;
    for (int i=0; i!=3; ++i){
      order.push_back(symMat[i][i]); 
      //std::cout << "   " << symMat[i][i] << std::endl;
      }
    sort(order.begin(), order.end());
    //double lambdaOne   = order[2]; 
    double lambdaTwo   = order[1];
    double lambdaThree = order[0];

    //std::cout << "sum of lambdas =  " << lambdaOne + lambdaTwo + lambdaThree << std::endl;

    sphericity_ = 3/2 * (lambdaTwo + lambdaThree); 
    aplanarity_ = 3/2 *  lambdaThree;
    planarity_  = 2 * (sphericity_ - 2*aplanarity_) /3; 

  }
  else { cerr << "Error: momentum tensor not symmetric" << std::endl; }

}
