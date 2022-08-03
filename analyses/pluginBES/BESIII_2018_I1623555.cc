// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  eta' -> eta pi0 pi0 or eta pi+ pi-
  class BESIII_2018_I1623555 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1623555);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==PID::ETAPRIME);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      ETA.addStable(PID::PI0);
      ETA.addStable(PID::K0S);
      ETA.addStable(PID::ETA);
      declare(ETA, "ETA");
      // histograms
      for(unsigned int ix=0;ix<4;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  if(ix<2)
	    book(_h[ix][iy],1+ix,1,1+iy);
	  else
	    book(_h[ix][iy],"TMP/h_"+toString(ix+1)+"_"+toString(iy+1),refData(1+ix,1,1+iy));
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { {221,1}, { 111,2} };
      static const map<PdgId,unsigned int> & mode2 = { {221,1}, { 211,1}, {-211,1} };
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for(unsigned int ix=0;ix<ETA.decaying().size();++ix) {
	// eta' -> 2pi0 eta
	if ( ETA.modeMatches(ix,3,mode1) ) {
	  const Particle &  eta = ETA.decayProducts()[ix].at(221)[0];
	  const Particles & pi0 = ETA.decayProducts()[ix].at(111);
	  double s1 = (eta   .momentum()+pi0[0].momentum()).mass2();
	  double s2 = (eta   .momentum()+pi0[1].momentum()).mass2();
	  double s3 = (pi0[0].momentum()+pi0[1].momentum()).mass2();
	  double mOut = eta.mass()+pi0[0].mass()+pi0[1].mass();
	  double Q = ETA.decaying()[ix].mass()-mOut;
	  double X = sqrt(3.)/2./ETA.decaying()[ix].mass()/Q*fabs(s1-s2);
	  double Y = mOut/2./Q/pi0[0].mass()/ETA.decaying()[ix].mass()*
	    (sqr(ETA.decaying()[ix].mass()-eta.mass())-s3)-1.;
	  _h[1][0]->fill(X);
	  _h[1][1]->fill(Y);
	  _h[2][0]->fill(sqrt(s3));
	  _h[3][0]->fill(sqrt(s1));
	  _h[3][0]->fill(sqrt(s2));
	}
	// eta' -> pi+ pi- eta
	else if ( ETA.modeMatches(ix,3,mode2) ) {
	  const Particle & eta = ETA.decayProducts()[ix].at( 221)[0];
	  const Particle & pip = ETA.decayProducts()[ix].at( 211)[0];
	  const Particle & pim = ETA.decayProducts()[ix].at(-211)[0];
	  double s1 = (eta.momentum()+pim.momentum()).mass2();
	  double s2 = (eta.momentum()+pip.momentum()).mass2();
	  double s3 = (pim.momentum()+pip.momentum()).mass2();
	  double mOut = eta.mass()+pip.mass()+pim.mass();
	  double Q = ETA.decaying()[ix].mass()-mOut;
	  double X = sqrt(3.)/2./ETA.decaying()[ix].mass()/Q*(s1-s2);
	  double Y = mOut/2./Q/pip.mass()/ETA.decaying()[ix].mass()*
	    (sqr(ETA.decaying()[ix].mass()-eta.mass())-s3)-1.;
	  _h[0][0]->fill(X);
	  _h[0][1]->fill(Y);
	  _h[2][1]->fill(sqrt(s3));
	  _h[3][1]->fill(sqrt(s1));
	  _h[3][1]->fill(sqrt(s2));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      vector<double> phsp[2][2] = {{{0.4958173595002211,1.9619294234849218,2.769343347436134,3.3754760296363586,3.878496706355459,4.314870895705234,4.70312212702334,
				     5.054215323586935,5.375279853746009,5.671257584168677,5.9457370820689865,6.201417631268466,6.44038618211146,6.664291971724481,
				     6.874461576204892,7.071977500668689,7.257733504274686,7.432474552616055,7.596826304190006,7.751317284426186,7.896395833035851,
				     8.032443239272407,8.15978404588282,8.278694215122778,8.389407655537537,8.492121473687217,8.587000220318588,8.674179332711457,
				     8.753767925627045,8.825851046861633,8.89049148607537,8.947731204719013,8.997592438695575,9.040078512568408,9.075174393713954,
				     9.102847006109055,9.123045315893089,9.135700194014767,9.140724054787638,9.13801026269368,9.127432292951557,9.108842623839646,
				     9.082071330112504,9.046924336557808,9.003181278149844,8.950592897494763,8.888877890152871,8.81771908233236,8.73675879108212,
				     8.645593171181584,8.543765290652805,8.430756591127544,8.305976269528442,8.168747947286356,8.018292746910289,7.853707531862115,
				     7.673936516314864,7.477733601051388,7.26361143821108,7.029771003081123,6.774001658321787,6.493534946280308,6.184822709548511,
				     5.843185002115959,5.462219359013247,5.03273598213602,4.540644424927434,3.962141200519649,3.2500824294265236,2.274198460941213},
				    {1.3100886831082426,2.59729040074684,3.381299983575231,4.000787719010948,4.524937634762839,4.983978373623454,5.3944278031342705,
				     5.7664948628407044,6.107010795600696,6.4208084604346265,6.711450584339649,6.981647651964704,7.233513384923498,7.468728930123908,
				     7.688652731356644,7.894396548261923,8.08687954469034,8.266867695860503,8.435003085521664,8.591826066910466,8.737792275150909,
				     8.87328585137746,8.998629828980688,9.114094358145204,9.219903257387287,9.31623925009379,9.403248151226878,9.481042202232233,
				     9.549702702777214,9.609282050887044,9.65980527467064,9.701271116566254,9.733652713050626,9.756897897684288,9.770929142161672,
				     9.775643137844083,9.770910008327416,9.756572131216842,9.732442533692973,9.698302810793034,9.653900496532554,9.598945794679436,
				     9.533107546305716,9.456008272606084,9.367218080231398,9.266247147243595,9.152536412939227,9.025445962449908,8.884240409151237,
				     8.728070306173185,8.555948217222312,8.366717471550364,8.15901069111041,7.931193686998973,7.681287872309445,7.4068601618098855,
				     7.10486189676173,6.771384421838216,6.40127126802985,5.987467571923839,5.519847572730394,4.982887990977129,4.3503709370253265,
				     3.5703872094268094,2.4995596481282525}},
				   {{0.629013798471671,2.1866118268994468,3.0368605864518217,3.6771490374606355,4.207007718753613,4.66452836989242,5.069326043954977,
				     5.433115314702656,5.763565023040947,6.066021922820969,6.344389357111023,6.601618652129766,6.8400036740352945,7.061367217314588,
				     7.26718418933565,7.458665984819544,7.636820030022579,7.8024928814423395,7.95640210521801,8.099160303672743,8.231293520219943,
				     8.353255538761486,8.46543913064336,8.568184994935699,8.661788929350681,8.74650762487693,8.82256337558689,8.890147922266685,
				     8.949425595529279,9.000535884920497,9.043595531165272,9.078700216339863,9.105925909443517,9.12532991118698,9.136951630805799,
				     9.14081311859039,9.13691937002785,9.125258410513581,9.105801163130566,9.078501095666466,9.043293636510745,9.000095341996824,
				     8.94880278973171,8.889291163010668,8.821412479942126,8.744993406613599,8.659832575438015,8.565697306248044,8.462319596663876,
				     8.34939120675572,8.226557606697057,8.093410478498992,7.949478354280492,7.794214818859097,7.626983480150643,7.4470385789742775,
				     7.253499607649825,7.045317527593014,6.821228932575126,6.5796924547567555,6.318798205008482,6.036134778147875,5.728586585466067,
				     5.392010752627947,5.020692096532121,4.606354301328156,4.136181095725286,3.5882580110647293,2.918423656686351,1.9999831254556695},
				    { 0.23873578916690366,2.0500489551803267,3.1201857892092946,3.8789686604099263,4.493437017697335,5.017527030470675,5.477296385866756,
				      5.887761534500594,6.258527506968923,6.596182009973793,6.905475510266001,7.189965315739698,7.452394570462295,7.694928593219781,
				      7.919309232207209,8.126959552432252,8.319057099733737,8.496586547202293,8.660378387673369,8.811137925405763,8.949467363242505,
				      9.075882871818438,9.190827942723814,9.294683942100447,9.387778521184131,9.470392361301753,9.542764605201127,9.605097236829916,
				      9.657558606430447,9.70028624956067,9.733389112310816,9.756949267063849,9.771023181187681,9.775642583275978,9.770814956656556,
				      9.756523676842061,9.732727797572638,9.699361478358194,9.65633303426016,9.603523575310243,9.540785187582555,9.467938589449417,
				      9.384770173568707,9.29102831578889,9.186418793789159,9.070599107145938,8.943171421166314,8.803673761129039,8.651568949209869,
				      8.48623058428439,8.306925084850821,8.11278839862622,7.902795347561311,7.6757185844345175,7.4300725389000055,7.16403506942705,
				      6.875334930613836,6.561084820602421,6.217523820876327,5.839600487397429,5.420255818696847,4.949088517083736,4.409586318075515,
				      3.772390577191644,2.973878032379262,1.7830364647510994,0.05380912975709777}}};

      // normalize histograms
      for(unsigned int ix=0;ix<4;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  normalize(_h[ix][iy]);
	}
      }
      // last two plots  convert to scatter and normalize to phase space volume in bin
      double step = 0.002;
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  // convert to scatter
	  Scatter2DPtr tmp;
	  book(tmp,3+ix,1,1+iy);
	  barchart(_h[ix+2][iy],tmp);
	  // divide by phase space volume
	  for(unsigned int ip=0;ip<tmp->points().size();++ip) {
	    tmp->points()[ip].scaleY(1./phsp[ix][iy][ip]/step);
	  }
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1623555);

}
