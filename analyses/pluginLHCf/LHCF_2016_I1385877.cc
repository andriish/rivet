// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableParticles.hh"
namespace Rivet {


/// @brief Add a short analysis description here
class LHCF_2016_I1385877 : public Analysis {
public:

	/// Constructor
	DEFAULT_RIVET_ANALYSIS_CTOR(LHCF_2016_I1385877);

	//In case of some models there can be very small value pT but greater than 0.
	//In order to avoid unphysical behavior in the first bin a cutoff is needed
	//If you are sure the model does not have this problem you can set pt_cutoff to 0.
	const double pt_cutoff = 0.01;

	/// @name Analysis methods
	//@{
        void bookHistosPP(const double _sqrts) {
            if (fuzzyEquals( _sqrts/GeV, 7000., 1E-3)) {
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(  8.8,  9.0, book(tmp,  2, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(  9.0,  9.2, book(tmp,  3, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(  9.2,  9.4, book(tmp,  4, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(  9.4,  9.6, book(tmp,  5, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(  9.6,  9.8, book(tmp,  6, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(  9.8, 10.0, book(tmp,  7, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 10.0, 10.2, book(tmp,  8, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 10.2, 10.4, book(tmp,  9, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 10.4, 10.6, book(tmp, 10, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 10.6, 10.8, book(tmp, 11, 1, 2));}

                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.0, 0.2, book(tmp, 12, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.2, 0.4, book(tmp, 13, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.4, 0.6, book(tmp, 14, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.6, 0.8, book(tmp, 15, 1, 2));}
                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.8, 1.0, book(tmp, 16, 1, 2));}

                book(_p_pi0_rap_apT,      1, 1, 2);
                book(_h_pi0_rap,         21, 1, 2);
                book(_p_pi0_raploss_apT, 22, 1, 2);
                book(_h_pi0_raploss,     23, 1, 2);
            } 
            else if (fuzzyEquals( _sqrts/GeV, 2760., 1E-3)){
                    book(_p_pi0_rap_apT, 1, 1, 1);

                    {Histo1DPtr tmp; _h_pi0_rap_pT.add(  8.8, 9.0, book(tmp, 2, 1, 1));}
                    {Histo1DPtr tmp; _h_pi0_rap_pT.add(  9.0, 9.2, book(tmp, 3, 1, 1));}
                    {Histo1DPtr tmp; _h_pi0_rap_pT.add(  9.2, 9.4, book(tmp, 4, 1, 1));}
                    {Histo1DPtr tmp; _h_pi0_rap_pT.add(  9.4, 9.6, book(tmp, 5, 1, 1));}
                    {Histo1DPtr tmp; _h_pi0_rap_pT.add(  9.6, 9.8, book(tmp, 6, 1, 1));}

                    {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.0, 0.2, book(tmp, 12, 1, 1));}
                    {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.2, 0.4, book(tmp, 13, 1, 1));}

                    book(_h_pi0_rap, 21, 1, 1);

                    book(_p_pi0_raploss_apT, 22, 1, 1);
                    book(_h_pi0_raploss, 23, 1, 1);
            }
            else{
                    MSG_WARNING("p-p collisions : energy out of range!");
            }
        }

        void bookHistosPPb(const double _sqrts) {
            if (fuzzyEquals(_sqrts/sqrt(208.)/GeV, 5020., 1E-3)) {

                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 8.8,  9.0, book(tmp,  2, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 9.0,  9.2, book(tmp,  3, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 9.2,  9.4, book(tmp,  4, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 9.4,  9.6, book(tmp,  5, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 9.6,  9.8, book(tmp,  6, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add( 9.8, 10.0, book(tmp,  7, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(10.0, 10.2, book(tmp,  8, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(10.2, 10.4, book(tmp,  9, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(10.4, 10.6, book(tmp, 10, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_rap_pT.add(10.6, 10.8, book(tmp, 11, 1, 3));}

                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.0,  0.2, book(tmp, 12, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.2,  0.4, book(tmp, 13, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.4,  0.6, book(tmp, 14, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.6,  0.8, book(tmp, 15, 1, 3));}
                {Histo1DPtr tmp; _h_pi0_pT_pZ.add(  0.8,  1.0, book(tmp, 16, 1, 3));}

                book(_p_pi0_rap_apT,      1, 1, 3);
                book(_p_pi0_raploss_apT, 22, 1, 3);
            }
            else {
                    MSG_WARNING("p-Pb collisions : energy out of range!");
            }
        }

	/// Book histograms and initialise projections before the run
	void init() {

            // Initialise and register projections
            declare(UnstableParticles(), "UFS");
            declare(Beam(), "Beam");

            // calculate beam rapidity
            const Particle bm1 = beams().first;
            const Particle bm2 = beams().second;
            _beam_rap_1 = bm1.rap();
            _beam_rap_2 = bm2.rap();
            MSG_INFO("Beam 1 : momentum " << bm1.pz() << " PID " << bm1.pid() << " rapidity " << bm1.rap() );
            MSG_INFO("Beam 2 : momentum " << bm2.pz() << " PID " << bm2.pid() << " rapidity " << bm2.rap() );

            _beam_rap = _beam_rap_1;
            
            const double _sqrts = sqrtS();
            MSG_INFO("CM energy: " << _sqrts );
            // p-p
            if (bm1.pid()==2212 && bm2.pid()==2212) {
                _isPP = true;
                bookHistosPP(_sqrts);
            }
            // p-Pb
            else if (bm1.pid()==PID::PROTON && bm2.pid()==PID::LEAD) {
                _isPP = false;
                bookHistosPPb(_sqrts);
            } 
            else MSG_WARNING("Beam PDGID out of range --- should be pp or p-Pb");
            _nevt = 0.;
	}

	/// Perform the per-event analysis
	void analyze(const Event& event) {
            _nevt++;

            // Select neutral pions
            const UnstableParticles &ufs = applyProjection<UnstableParticles> (event, "UFS");
            const auto pions = ufs.particles(Cuts::pz>0 && Cuts::abspid == 111 && Cuts::pT > pt_cutoff*GeV);
            for (auto p: pions ) {
                const double pT = p.pT()/GeV;
                const double rap = p.rap();
                const double raploss = _beam_rap_1 - p.rap();
                _p_pi0_rap_apT    ->fill(    rap,  p.pT()/MeV,          1.0);
                _p_pi0_raploss_apT->fill(raploss,  p.pT()/MeV,          1.0);
                _h_pi0_rap_pT      .fill(    rap,          pT,       1.0/pT);
                _h_pi0_pT_pZ.       fill(     pT,  p.pz()/GeV, p.E()/GeV/pT);
                if (_isPP) {
                  _h_pi0_rap        ->fill(rap    , 1.0);
                  _h_pi0_raploss    ->fill(raploss, 1.0);
                }
            }
	}

	/// Normalise histograms etc., after the run
	void finalize() {
		const double inv_scale_factor = 1. / _nevt / (2.*PI);
		const double pt_bin_width = 0.2;

		for (Histo1DPtr h: _h_pi0_pT_pZ.histos()){
			if(h->path() == "/LHCF_2016_I1385877/d12-x01-y01" ||
					h->path() == "/LHCF_2016_I1385877/d12-x01-y02" ||
					h->path() == "/LHCF_2016_I1385877/d12-x01-y03")
				h->scaleW( inv_scale_factor / (pt_bin_width-pt_cutoff) );
			else
				h->scaleW( inv_scale_factor / pt_bin_width );
		}

		const double scale_factor =  1. / _nevt / (2.*PI);
		const double rap_bin_width = 0.2;
		for (Histo1DPtr h: _h_pi0_rap_pT.histos()) {
			const int cutoff_bin = h->binIndexAt(pt_cutoff);

			if(cutoff_bin>=0) {
				const double cutoff_wdt = h->bin(cutoff_bin).xMax()-h->bin(cutoff_bin).xMin();
				h->bin(cutoff_bin).scaleW((cutoff_wdt)/(cutoff_wdt-pt_cutoff));
			}

			h->scaleW( scale_factor / rap_bin_width );
		}

		if(_isPP) {
			scale( _h_pi0_rap , 1. / _nevt );
			scale( _h_pi0_raploss , 1. / _nevt );
		}
	}
	//@}


private:

	/// @name Histograms
	//@{
	bool _isPP;
	double _nevt;
	double _beam_rap;
	double _beam_rap_1;
	double _beam_rap_2;
	BinnedHistogram _h_pi0_pT_pZ;
	BinnedHistogram _h_pi0_rap_pT;
	Profile1DPtr _p_pi0_rap_apT;
	Histo1DPtr _h_pi0_rap;
	Profile1DPtr _p_pi0_raploss_apT;
	Histo1DPtr _h_pi0_raploss;
	//@}
};


// The hook for the plugin system
DECLARE_RIVET_PLUGIN(LHCF_2016_I1385877);

}
