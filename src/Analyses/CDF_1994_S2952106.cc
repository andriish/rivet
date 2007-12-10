// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_1994_S2952106.hh"
using namespace Rivet;

#include "Rivet/RivetAIDA.hh"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////


// Book histograms
void CDF_1994_S2952106::init() {
  // Use histogram auto-booking
  //_histJetAzimuth_pTmax75_100  = bookHistogram1D(1, 2, 1, "Jet Jet azimuthal angle, pTmax=75..100");
  //_histJetAzimuth_pTmax100_130 = bookHistogram1D(2, 2, 1, "Jet Jet azimuthal angle, pTmax=100..130");
  //_histJetAzimuth_pTmax130_180 = bookHistogram1D(3, 2, 1, "Jet Jet azimuthal angle, pTmax=130..180");
  //_histJetAzimuth_pTmax180_    = bookHistogram1D(4, 2, 1, "Jet Jet azimuthal angle, pTmax>180");

  //const string hname = "HvsDphi";
  //const string htitle = "H vs Delta phi";
  //_histHvsDphi = bookHistogram2D(hname, htitle, 40, -4., 4., 32, 0., 3.2);

  //const string hname2 = "RvsAlpha";
  //const string htitle2 = "R vs alpha";
  //_histRvsAlpha = bookHistogram2D(hname2, htitle2, 50, 0., 5., 32, -1.6, 1.6);

  const string hname3 = "Jet1Et";
  const string htitle3 = "Et of leading jet";
  _histJet1Et = bookHistogram1D(hname3, htitle3, 40, 0., 500.);
  
  const string hname4 = "Jet2Et";
  const string htitle4 = "Et of 2nd leading jet";
  _histJet2Et = bookHistogram1D(hname4, htitle4, 40, 0., 500.);
  
  const string hname5 = "R23";
  const string htitle5 = "R distance between 2nd, 3rd jet";
  _histR23 = bookHistogram1D(hname5, htitle5, 50, 0., 5.);
  
  const string hname6 = "Jet3eta";
  const string htitle6 = "Pseudorapidity of 3rd jet";
  _histJet3eta = bookHistogram1D(hname6, htitle6, 42, -4., 4.);
  
  const string hname7 = "alpha";
  const string htitle7 = "alpha";
  _histAlpha = bookHistogram1D(hname7, htitle7, 42, -PI/2., PI/2.);
  
  //const string hname8 = "alphaMCvsDat";
  //const string htitle8 = "alpha MC vs. Data ";
  //_histAlphaMCvsDat = bookHistogram1D(hname8, htitle8, 42, -PI/2., PI/2.);


  const string hname9 = "alphaIdeal";
  const string htitle9 = "alpha Ideal";
  _histAlpaIdeal = bookHistogram1D(hname9, htitle9, 42, -PI/2., PI/2.);

  const string hname10 = "alphaCDF";
  const string htitle10 = "alpha CDF";  
  _histAlpaCDF = bookHistogram1D(hname10, htitle10, 42, -PI/2., PI/2.);

  const string hname11 = "R23Ideal";
  const string htitle11 = "R23 Ideal";  
  _histR23Ideal = bookHistogram1D(hname11, htitle11, 50, 0., 5.);

  const string hname12 = "R23CDF";
  const string htitle12 = "R23 CDF";  
  _histR23CDF = bookHistogram1D(hname12, htitle12, 50, 0., 5.);

  const string hname13 = "Jet3etaIdeal";
  const string htitle13 = "Jet3 eta Ideal";  
  _histJet3etaIdeal = bookHistogram1D(hname13, htitle13, 42, -4., 4.);

  const string hname14 = "Jet3etaCDF";
  const string htitle14 = "Jet3 eta CDF";  
  _histJet3etaCDF = bookHistogram1D(hname14, htitle14, 42, -4., 4.);




  _eventsTried = 0.0;

}


// Do the analysis
void CDF_1994_S2952106::analyze(const Event & event) {
  Log& log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;

  //increment counter for the number of events analysed
  _eventsTried++;
  //cout << "try event # " << int(_eventsTried) << endl;

  // Analyse and print some info  
  #ifdef HAVE_FASTJET
  const FastJets& jetpro = event.applyProjection(_conejetsproj);
  #else
  const D0ILConeJets& jetpro = event.applyProjection(_conejetsproj);
  #endif

  log << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.getNJets() << endl;

  // Find vertex and check  that its z-component is < 60 cm from the nominal IP
  const PVertex& pv = event.applyProjection(_vertexproj);
  /// @todo z- value assumed to be in mm, PYTHIA convention: dangerous!
  if (fabs(pv.getPrimaryVertex().position().z()) < 600.0) {

    const TotalVisibleMomentum& caloMissEt = event.applyProjection(_calmetproj);
    log << Log::DEBUG << "CaloMissEt.getMomentum().pT() = " << caloMissEt.getMomentum().pT() << endl;
    if (caloMissEt.getMomentum().pT()/sqrt(caloMissEt.getSET()) < 6.0) {

      /// @todo Remove compile-time flags like this when possible.
      #ifdef HAVE_FASTJET
      typedef vector<fastjet::PseudoJet> Jets;
      const Jets& jets = jetpro.getJetsPt();
      #else
      typedef list<FourMomentum> Jets;
      const Jets& jets = jetpro.getLorentzJets();
      #endif
 
      list<FourMomentum>::const_iterator jet1stPt = jets.end();
      list<FourMomentum>::const_iterator jet2ndPt = jets.end();
      list<FourMomentum>::const_iterator jet3rdPt = jets.end();
      log << Log::DEBUG << "jetlist size = " << jets.size() << endl;
      
      int Njet = 0;
      int NjetPtCut = 0;
      for (list<FourMomentum>::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
	log << Log::DEBUG << "List item pT = " << jt->pT() << " E=" << jt->E() 
	    << " pz=" << jt->pz() << endl;
	if (jt->pT() > 100.0) ++NjetPtCut;
	log << Log::DEBUG << "Jet pT =" << jt->pT() << " y=" << jt->pseudorapidity() 
	    << " phi=" << jt->azimuthalAngle() << endl; 
	if (jet1stPt == jets.end() || jt->pT() > jet1stPt->pT()) {
	  jet3rdPt = jet2ndPt;
	  jet2ndPt = jet1stPt;
	  jet1stPt = jt;
	  Njet++;
	} else if (jet2ndPt == jets.end() || jt->pT() > jet2ndPt->pT()) {
	  jet3rdPt = jet2ndPt;
	  jet2ndPt = jt;
	  Njet++;
	} else if (jet3rdPt == jets.end() || jt->pT() > jet3rdPt->pT()) {
	  jet3rdPt = jt;
	  Njet++;
	}
      }
      
      if (NjetPtCut >= 1) {
	log << Log::DEBUG << "Jet multiplicity after pT > 100GeV cut = " << NjetPtCut << endl; 
      }
      
      if (Njet>=3 && fabs(jet1stPt->pseudorapidity())<0.7 && fabs(jet2ndPt->pseudorapidity())<0.7) {
	log << Log::DEBUG << "Jet eta and pT requirements fulfilled" << endl;
	
	if (fabs(fabs(jet1stPt->azimuthalAngle()-jet2ndPt->azimuthalAngle())-PI) < 1./9.*PI/2.) {
	  log << Log::DEBUG << "1st & 2nd Jet phi requirement fulfilled" << endl;
	  
	  
	  _histJet1Et->fill(jet1stPt->pT(), event.weight());
	  _histJet2Et->fill(jet2ndPt->pT(), event.weight());
	  //_histR23->fill(jet2ndPt->deltaR(*jet3rdPt), event.weight());
	  _histR23->fill(delta_rad(jet2ndPt->pseudorapidity(),jet2ndPt->azimuthalAngle(), 
			      jet3rdPt->pseudorapidity(),jet3rdPt->azimuthalAngle() ), 
			 event.weight());
	  _histJet3eta->fill(jet3rdPt->pseudorapidity(), event.weight());
	  
			 
	  //next cut only required for alpha studies
	  if (jet3rdPt->pT() > 10.) {
	    log << Log::DEBUG << "jet3rdPt->pT()=" << jet3rdPt->pT() 
		 << " (>10.)" << endl;

	    double dPhi = fabs(jet3rdPt->azimuthalAngle() - jet2ndPt->azimuthalAngle());
	    dPhi -= int(dPhi/PI); //dPhi % PI (modulo)
	    
	    double dH = jet3rdPt->pseudorapidity() - jet2ndPt->pseudorapidity();
	    if (jet2ndPt->pseudorapidity() < 0.) dH *= -1;
	    
	    double alpha = atan(dH/dPhi);
	    //cout << "alpha=" << alpha << "  dH=" << dH 
	    // << "  dPhi=" << dPhi << endl;

	    _histAlpha->fill(alpha, event.weight());
	    //_histAlphaMCvsDat->fill(alpha, event.weight());

	  } //3rd jet Pt cut
        } //delta phi 1st, 2nd jet
      } //1st + 2nd jet pseudoRapidity & Njet>=3
    } //MET/sqrt(SET) cut
  } //z-vertex
  
  
  // Finished
  log << Log::DEBUG << "Finished analyzing" << endl;
}


// Finalize
void CDF_1994_S2952106::finalize() { 
  Log& log = getLog();

  /*
  /// @todo Apply correction
  double a, b, c, erra, errb, errc;
  for (int ibin = 0;  ibin < _histAlpha->getNbins(); ++ibin) {
    a = _histAlpha->GetBinContent(ibin);
    erra = _histAlpha->GetBinError(ibin);
    b = _histAlpaIdeal->GetBinContent(ibin);
    errb = _histAlpaIdeal->GetBinError(ibin);
    c = _histAlpaCDF->GetBinContent(ibin);
    errc = _histAlpaCDF->GetBinError(ibin);
    _histAlpha->SetBinContent(ibin, b/c);
    _histAlpha->SetBinError(ibin, sqrt(sqr(b)/sqr(c)*sqr(erra) + sqr(a)/sqr(c)*sqr(errb) + 
				       sqr(a*b/(sqr(c)))*sqr(errc) ) );
  }
  /// @todo Same correction to be applied for _hisR23 and _histJet3eta histograms
  */


  //normalise histograms to integrated publication luminosity of 4.2pb^-1    
  //crossSection function returns pb.
  double fac = 4.2 * crossSection() / _eventsTried;
  
  fac = 1.; //temp, until it works

  log << Log::INFO << "crossSection()=" << crossSection()
      << "  _eventsTried=" << _eventsTried 
      << "  scale factor=" << fac << endl;

  
  //AIDA::IHistogram2D* _histHvsDphi;
  //AIDA::IHistogram2D* _histRvsAlpha;
  _histJet1Et->scale(fac);
  _histJet2Et->scale(fac);
  _histR23->scale(fac);
  _histJet3eta->scale(fac);
  _histAlpha->scale(fac);
  //AIDA::IHistogram1D* _histAlphaMCvsDat;
  
}
