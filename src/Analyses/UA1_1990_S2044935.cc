#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Analyses/UA1_1990_S2044935.hh"

namespace Rivet {


    /// Default constructor
    UA1_1990_S2044935::UA1_1990_S2044935() 
    : Analysis("UA1_1990_S2044935")
   {
      setBeams(PROTON, ANTIPROTON);
      const ChargedFinalState cfs(-2.5, 2.5);
      const FinalState fs2(-6., 6.);
      const FinalState fs(-2.5,2.5);
      addProjection(fs, "FS");
      addProjection(fs2, "FS2");
      addProjection(ChargedFinalState(-2.5, 2.5), "CFS");
      addProjection(Beam(), "Beam");
      addProjection(TotalVisibleMomentum(fs), "Mom");
    }

    /// @name Analysis methods
    //@{
    //book histograms
    void UA1_1990_S2044935::init() { 

      _hist_sigma200 =
        bookHistogram1D(1,1,1);
      _hist_sigma500 =
        bookHistogram1D(1,1,2);
      _hist_sigma900 =
        bookHistogram1D(1,1,3);
      _hist_Esigma200 =
        bookHistogram1D(2,1,1);
      _hist_Esigma500 =
        bookHistogram1D(2,1,2);
      _hist_Esigma900 =
        bookHistogram1D(2,1,3);
      _hist_Esigmapoint8 =
        bookHistogram1D(3,1,1);
      _hist_Esigma4 =
        bookHistogram1D(4,1,1);
      _hist_Esigma8 =
        bookHistogram1D(5,1,1);
      _hist_Et200 =
        bookHistogram1D(9,1,1);
      _hist_Et500 =
        bookHistogram1D(10,1,1);
      _hist_Et900 =
        bookHistogram1D(11,1,1);
      _hist_Pt63 =
        bookProfile1D(8,1,1);
      _hist_Pt200 =
        bookProfile1D(6,1,1);
      _hist_Pt900 =
        bookProfile1D(7,1,1);
      _hist_Etavg200 =
        bookProfile1D(12,1,1); 
      _hist_Etavg500 =
        bookProfile1D(12,1,2); 
      _hist_Etavg900 =
        bookProfile1D(12,1,3);
    }
    
    void UA1_1990_S2044935::analyze(const Event& event) 
    {
      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
      const double weight = event.weight();
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      double multi = cfs.particles().size();

	if (fuzzyEquals(sqrtS/GeV, 200, 1E-4)) 
	{
	  _hist_sigma200->fill(multi, weight);
	} 
	else if (fuzzyEquals(sqrtS/GeV, 500)) 
	{
	  _hist_sigma500->fill(multi, weight);
	}
	else if (fuzzyEquals(sqrtS/GeV, 900)) 
	{
	  _hist_sigma900->fill(multi, weight);
	}
      foreach (const Particle& p, fs.particles())
      {
///@todo figure out where the extra factor of 0.5 comes from in the weight factor (eta range?).
        double pt = p.momentum().pT();
        if (fuzzyEquals(sqrtS/GeV, 200, 1E-4))
        {
          _hist_Esigma200->fill(pt, weight/(2.*10.*M_PI*pt));
        }
        if (fuzzyEquals(sqrtS/GeV, 500))
        {
          _hist_Esigma500->fill(pt, weight/(2.*10.*M_PI*pt));
        }
        if (fuzzyEquals(sqrtS/GeV, 900))
        {
          _hist_Esigma900->fill(pt, weight/(2.*10.*M_PI*pt));
          const double dnch_deta = multi/5.0;
          if (dnch_deta >= 0.8 && dnch_deta <= 4)
          {
            _hist_Esigmapoint8->fill(pt, weight/(10.*M_PI*pt));
          }
          else if (dnch_deta > 4 && dnch_deta <= 8)
          {
            _hist_Esigma4->fill(pt, weight/(10.*M_PI*pt));
          }
          else if(dnch_deta > 8)
          {
            _hist_Esigma8->fill(pt, weight/(10.*M_PI*pt));
          }
        }                
      }

      const double Et = applyProjection<TotalVisibleMomentum>(event, "Mom").scalarET();

	if (fuzzyEquals(sqrtS, 200, 1E-4)) 
	{
	  _hist_Et200->fill(Et, weight);
	} 
	else if (fuzzyEquals(sqrtS, 500)) 
	{
	  _hist_Et500->fill(Et, weight);
	}
	else if (fuzzyEquals(sqrtS, 900)) 
	{
	  _hist_Et900->fill(Et, weight);
	}

      foreach (const Particle& p, cfs.particles())
      {
        {
      	if (fuzzyEquals(sqrtS, 63, 1E-3)) 
	{
	  _hist_Pt63->fill(multi, p.momentum().pT(), weight);
	} 
	else if (fuzzyEquals(sqrtS, 200, 1E-4)) 
	{
	  _hist_Pt200->fill(multi, p.momentum().pT(), weight);
	  _hist_Etavg200->fill(multi, Et, weight);
	}
	else if (fuzzyEquals(sqrtS, 500))
	{
	  _hist_Etavg500->fill(multi, Et, weight);
	}
	else if (fuzzyEquals(sqrtS, 900)) 
	{
	  _hist_Pt900->fill(multi, p.momentum().pT(), weight);
	  _hist_Etavg900->fill(multi, Et, weight);
	}
      }
    }
    }


    
    void UA1_1990_S2044935::finalize() {
    ///@todo: get the total cross-sections from the generator
    ///@todo: check if the scaling for Esigmpoint8, Esigma4 and Esigma8 are correct.
       normalize(_hist_sigma200, 27.9);
      normalize(_hist_sigma500, 31.5);
      normalize(_hist_sigma900, 34.4);
       scale(_hist_Esigma200, 27.9/sumOfWeights());
     scale(_hist_Esigma500, 31.5/sumOfWeights());
      scale(_hist_Esigma900, 34.4/sumOfWeights());
      scale(_hist_Esigmapoint8, 34400./sumOfWeights());
      scale(_hist_Esigma4, 3440./sumOfWeights());
      scale(_hist_Esigma8, 344./sumOfWeights());
            normalize(_hist_Et200, 27.9);
      normalize(_hist_Et500, 31.5);
      normalize(_hist_Et900, 34.4);
    }



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<UA1_1990_S2044935> plugin_UA1_1990_S2044935;

}
