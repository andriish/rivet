// -*- C++ -*-

#include<iostream>

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/HepEx0107012.hh" 
using namespace Rivet;

#include "Rivet/RivetAIDA.hh"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


#define PI (3.1415926)




// Booking the histograms:
void HepEx0107012::init() {
	
  Log& log = getLog();
  log << Log::INFO << "HepEx0107012::init(): fired" << endl;

  
//ROOT histograms:  
  _h_pt_e    = bookHistogram1D("pt_e", "pt_e", 20, 0., 200.);
  //_h_pt_e    = bookHistogram1D(20, 1, 1, "pt_e");

  _h_pt_miss = bookHistogram1D("pt_miss", "pt_miss", 20, 0., 200.);
  //_h_pt_miss = bookHistogram1D(20, 1, 1, "pt_miss");

  string h_pt_w_title = "pt_W";
  //_h_pt_w = bookHistogram1D(23, 1, 1, h_pt_w_title);
  _h_pt_w = bookHistogram1D(h_pt_w_title, h_pt_w_title, _bins_pt_w);
  //_h_pt_w = bookHistogram1D(h_pt_w_title, h_pt_w_title, 20, 0., 200.);
  
  string h_pt_w_true_title = "pt_W_true";
  //_h_pt_w_true = bookHistogram1D(23, 1, 1, h_pt_w_true_title);
  _h_pt_w_true = bookHistogram1D(h_pt_w_true_title, h_pt_w_true_title, _bins_pt_w);
  //_h_pt_w_true = bookHistogram1D(h_pt_w_true_title, h_pt_w_true_title, 20, 0., 200.);
  
  string h_w_rec_eff_title = "W_rec_eff";
  //_h_w_rec_eff_factor = bookHistogram1D(23, 1, 1, h_w_rec_eff_factor_title);
  _h_w_rec_eff_factor = bookHistogram1D(h_w_rec_eff_title, h_w_rec_eff_title, _bins_pt_w);
  //_h_w_rec_eff_factor = bookHistogram1D(h_w_rec_eff_title, h_w_rec_eff_title, 20, 0., 200.);
  

  log << Log::DEBUG << " init(): return" << endl;
}


// Event analysis:

void HepEx0107012::analyze(const Event & event) {
	Log& log = getLog();
	log << Log::DEBUG<< "Starting analyzing" << endl;
    
	const double weight = event.weight();

	for (GenEvent::particle_const_iterator pi = event.genEvent().particles_begin();
             pi != event.genEvent().particles_end(); 
	     ++pi) {
		HepPDT::ParticleID pInfo = (*pi)->pdg_id();
		if (pInfo.abspid() != 24)
			continue;
		HepMC::FourVector fv = (*pi)->momentum();
            	LorentzVector  pw_true(fv.px(), fv.py(), fv.pz(), fv.e());
		_h_pt_w_true->fill(pw_true.perp(),weight);
	}
	
	const FinalState& fs      = event.applyProjection(_fsproj);

	ParticleVector _theIsolChargedLeptons;

	LorentzVector p_tot(0.,0.,0.,0.);

	for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
		HepPDT::ParticleID pInfo = p->getPdgId();
            	HepMC::FourVector fv = p->getMomentum();
            	LorentzVector  pmom(fv.px(), fv.py(), fv.pz(), fv.e());
	    	if ( !( pInfo.abspid() == 12 || pInfo.abspid() == 14 || pInfo.abspid() == 16 ) ) 
			    p_tot = p_tot + pmom;
	    	if (pInfo.abspid() == 11) { // || pInfo.abspid() == 13) { // e+-, mu+-
			if ( pmom.perp() >= 25. /*25*/  && abs(pmom.eta()) <= 2.5   ) { // <=== cuts from HEPEX-0106027 (see HEPEX-0010026)
				double E_cone_0_4 = 0., E_em_02 = 0., f_iso;
			    	log << Log::DEBUG << "e+-: " << pInfo.pid() << " " << pmom 
					         << std::endl << "-----" << std::endl;
			    	for (ParticleVector::const_iterator p1 = fs.particles().begin(); p1 != fs.particles().end(); ++p1) {
					HepPDT::ParticleID p1Info = p1->getPdgId();
				    	HepMC::FourVector fv1 = p1->getMomentum();
				    	LorentzVector  p1mom(fv1.px(), fv1.py(), fv1.pz(), fv1.e());
				    	string thesame = " ";   string in04    = " ";  
					double dE = 0.;
				    	if ( p1 == p ) 
						thesame = " thesame";
				    	if ( pmom.deltaR(p1mom) <= 0.4 ) {  
						in04 = " <= 0.4 ";  
						dE = p1mom.e();  
					}
			            	log << Log::DEBUG  << "FS particle: " << p1Info.pid() 
							  << " "    << p1mom << " "
							  << "dR: " << pmom.deltaR(p1mom) 
							  << thesame 
							  << in04 
							  << dE ;
				    	if ( p1Info.abspid() == 12 || p1Info.abspid() == 14 || p1Info.abspid() == 16 ) {
						log << Log::DEBUG << std::endl;
						continue; // invisible
					}
				    	if ( pmom.deltaR(p1mom) > 0.4) {
						log << Log::DEBUG << std::endl;
						continue; // <== cut: off 0.4 cone
				    	}
				    	E_cone_0_4 += dE;
				    	log << Log::DEBUG<< " | " << E_cone_0_4 << std::endl;
					
			    	}
				//FIXME: should be (E_tot(0.4) - E_em(0.2)) / E_em(0.2),
				// we use a simplified definition instead:
			    	f_iso = (E_cone_0_4 - pmom.e()) / pmom.e(); 
			    	if (f_iso < 0.15) { // <== cut
					_theIsolChargedLeptons.push_back(Particle(*p)); 
				    	log << Log::DEBUG<< std::endl
				                          << "Isolated e+-: " << pmom 
							  << " pt,e: " << pmom.perp()  << ","  <<  pmom.e() 
				        		  << " E_cone_0_4: " << E_cone_0_4
				        	 	  << " f_iso: " << f_iso 
				        		  << std::endl;
			    	}
		    	}
	    	} // End if (... abspid ... == 12 ...
    	} // End of cycle over FS particles
    

	double PTW = 0.;
	double pt_miss = p_tot.perp();
    
   	 _h_pt_miss->fill (pt_miss, weight);
   /**/ 
    	if (pt_miss < 25.) // cut: PT_miss >= 25 GeV
		return;    
   /**/

    	int nlep = _theIsolChargedLeptons.size();
    	//if (nlep == 2) { //TODO: what if there are 2 l+-'s not consistent with Z?
    	if (nlep == 2 && _theIsolChargedLeptons[0].getPdgId() == -_theIsolChargedLeptons[1].getPdgId() ) { //two same flavour leptons, opposite charge
		LorentzVector psum(0.,0.,0.,0.);
		for (ParticleVector::const_iterator p = _theIsolChargedLeptons.begin(); p != _theIsolChargedLeptons.end(); ++p) {
			HepMC::FourVector fv = p->getMomentum();
			LorentzVector  plep(fv.px(), fv.py(), fv.pz(), fv.e());
			psum += plep;
		}
		if ((psum.m()-91.18) < 15.) // Z candidate
		  return; // drop the event
	}
	else if (nlep == 1) {
	        Particle& p = _theIsolChargedLeptons[0];
		HepMC::FourVector fv = p.getMomentum();
		LorentzVector plep(fv.px(), fv.py(), fv.pz(), fv.e());
		LorentzVector PTW (plep.px()-p_tot.px(), plep.py()-p_tot.py(),0.,0.);
		log << Log::DEBUG << " PTW: " << PTW.perp() << std::endl;
		_h_pt_e->fill(plep.perp(), weight);
		_h_pt_w->fill(PTW.perp(), weight);
	}
    

    log << Log::DEBUG << "Finished analyzing" << std::endl;
    return;
}



void HepEx0107012::finalize() { 
  Log& log = getLog();

  _h_w_rec_eff_factor = histogramFactory().divide("/HepEx0107012/W_rec_eff", *_h_pt_w_true, *_h_pt_w); // true/measured(after cuts) => factors large one!

  //histogramFactory().destroy(_h_pt_w);
  //histogramFactory().destroy(_h_pt_w_true);

}


