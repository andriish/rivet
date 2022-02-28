// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief strange fraction at 7,8,13 TeV
  class LHCB_2019_I1760257 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2019_I1760257);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      //histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_c_B[ix],"TMP/c_B_"+toString(ix+1));
	book(_h_pT[ix],"TMP/h_pT_"+toString(ix+1),refData(5,1,1));
	for(unsigned int iy=0;iy<3;++iy) {
	  book(_h_pL[ix][iy],"TMP/h_pL_"+toString(ix+1)+"_"+toString(iy+1),refData(3,1,1+iy));
	}
	for(unsigned int iy=0;iy<5;++iy) {
	  book(_h_kin[ix][iy],"TMP/h_kin_"+toString(ix+1)+"_"+toString(iy+1),refData(4,1,1+iy));
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for( const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==521 or
      										Cuts::abspid==531)) {
      	// kinematic region
      	// rapidity cuts
      	double y = p.rapidity(), eta = p.eta();
      	if(y<2. || y> 4.5 || eta<2. || eta>6.5) continue;
      	// momentum cuts
      	double pT = p.perp(), pL = p.momentum().z(), pTot=p.momentum().p3().mod();
      	if(pT<0.5 || pT>40. || pTot<20. || pTot>700. || pL<20. || pL>700.) continue;
      	// type of meson
      	unsigned int imeson = (p.abspid()%100)/10 - 2;
      	_c_B[imeson]->fill();
      	_h_pT[imeson]->fill(pT);
      	if (pL<75.)       _h_pL[imeson][0]->fill(pT);
      	else if (pL<125.) _h_pL[imeson][1]->fill(pT);
      	else if (pL>700.) _h_pL[imeson][2]->fill(pT);
      	_h_kin[imeson][0]->fill(pTot);
      	_h_kin[imeson][1]->fill(pL  );
      	_h_kin[imeson][2]->fill(pT  );
      	_h_kin[imeson][3]->fill(eta );
      	_h_kin[imeson][4]->fill(y   );
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // ratio of branching ratios from PDSG 2020
      double brRatio = 1.08e-3/1.02e-3;
      for(unsigned int iy=0;iy<5;++iy) {
      	Scatter2DPtr tmp;
      	book(tmp,4,1,iy+1);
      	divide(_h_kin[1][iy],_h_kin[0][iy],tmp);
	tmp->scaleY(brRatio);
      	if(iy>=3) continue;
      	book(tmp,3,1,iy+1);
      	divide(_h_pL [1][iy],_h_pL [0][iy],tmp);
	tmp->scaleY(brRatio);
      }
      unsigned int is=0;
      if (isCompatibleWithSqrtS(7000)) {
      	is=0;
      }
      else if (isCompatibleWithSqrtS(8000)) {
      	is=1;
      }
      else if (isCompatibleWithSqrtS(13000)) {
      	is=2;
      }
      else  {
      	throw Error("Invalid CMS energy for LHCB_2019_I1760257");
      }
      Scatter2DPtr tmp;
      book(tmp,5,1,1+is);
      divide(_h_pT[1],_h_pT[0],tmp);
      tmp->scaleY(brRatio);
      // ratio
      Scatter1D s1("/LHCB_2019_I1760257/TMP/r_1");
      Scatter1DPtr s1d = registerAO(s1);
      divide(_c_B[1],_c_B[0],s1d);
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr ratio;
      book(ratio, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
      	const double x  = temphisto.point(b).x();
      	pair<double,double> ex = temphisto.point(b).xErrs();
      	if( fuzzyEquals(x,sqrtS()/1000.,1e-2) )
      	  ratio->addPoint(x,s1d->points()[0].x(),ex,s1d->points()[0].xErrs());
      	else
      	  ratio->addPoint(x, 0., ex, make_pair(0.,.0));
      }
      ratio->scaleY(brRatio);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_B[2];
    Histo1DPtr _h_pL[2][3],_h_kin[2][5];
    Histo1DPtr _h_pT[2];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(LHCB_2019_I1760257);

}
