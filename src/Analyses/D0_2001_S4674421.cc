// -*- C++ -*-

#include<iostream>

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/D0_2001_S4674421.hh" 
using namespace Rivet;

#include "Rivet/RivetAIDA.hh"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


// Booking the histograms:
void D0_2001_S4674421::init() {	
  Log& log = getLog();
  log << Log::INFO << "D0_2001_S4674421::init(): fired" << endl;

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


void D0_2001_S4674421::analyze(const Event & event) {
	Log& log = getLog();
	log << Log::DEBUG<< "Starting analyzing" << endl;    
	const double weight = event.weight();
    
	for (GenEvent::particle_const_iterator pi = event.genEvent().particles_begin();
         pi != event.genEvent().particles_end(); ++pi) {
      HepPDT::ParticleID pInfo = (*pi)->pdg_id();
      if (pInfo.abspid() != 24) continue;
      FourMomentum pw_true((*pi)->momentum());
      _h_pt_w_true->fill(pw_true.polarRadius(), weight);
	}
	
	const FinalState& fs = event.applyProjection(_fsproj);
    
	ParticleVector theIsolChargedLeptons;
	FourMomentum p_tot;
	for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      HepPDT::ParticleID pInfo = p->getPdgId();
      FourMomentum pmom(p->getMomentum());
      if ( !( pInfo.abspid() == 12 || pInfo.abspid() == 14 || pInfo.abspid() == 16 ) ) p_tot = p_tot + pmom;
      if (pInfo.abspid() == 11) { // || pInfo.abspid() == 13) { // e+-, mu+-
        // Cuts from HEPEX-0106027 (see HEPEX-0010026)
        if ( pmom.polarRadius() >= 25. /*25*/  && fabs(pmom.pseudorapidity()) <= 2.5   ) {
          double E_cone_0_4 = 0., f_iso;
          //double E_em_02 = 0.;
          log << Log::DEBUG << "e+-: " << pInfo.pid() << " " << pmom 
              << endl << "-----" << endl;
          for (ParticleVector::const_iterator p1 = fs.particles().begin(); p1 != fs.particles().end(); ++p1) {
            HepPDT::ParticleID p1Info = p1->getPdgId();
            FourMomentum p1mom(p1->getMomentum());
            string thesame = " ";   string in04 = " ";  
            double dE = 0.;
            if ( p1 == p ) thesame = " thesame";
            if ( deltaR(pmom, p1mom) <= 0.4 ) {  
              in04 = " <= 0.4 ";
              dE = p1mom.E();
            }
            log << Log::DEBUG  << "FS particle: " << p1Info.pid() 
                << " "    << p1mom << " "
                << "dR: " << deltaR(pmom, p1mom)
                << thesame 
                << in04 
                << dE ;
            if ( p1Info.abspid() == 12 || p1Info.abspid() == 14 || p1Info.abspid() == 16 ) {
              log << Log::DEBUG << endl;
              continue; // invisible
            }
            if ( deltaR(pmom, p1mom) > 0.4) {
              log << Log::DEBUG << endl;
              continue; // <== cut: off 0.4 cone
            }
            E_cone_0_4 += dE;
            log << Log::DEBUG<< " | " << E_cone_0_4 << endl;
			
          }
          //FIXME: should be (E_tot(0.4) - E_em(0.2)) / E_em(0.2),
          // we use a simplified definition instead:
          f_iso = (E_cone_0_4 - pmom.E()) / pmom.E(); 
          if (f_iso < 0.15) { // <== cut
            theIsolChargedLeptons.push_back(Particle(*p)); 
            log << Log::DEBUG<< std::endl
                << "Isolated e+-: " << pmom 
                << " pt,e: " << pmom.polarRadius()  << ","  <<  pmom.E() 
                << " E_cone_0_4: " << E_cone_0_4
                << " f_iso: " << f_iso 
                << endl;
          }
        }
      } // End if (... abspid ... == 12 ...
    } // End of cycle over FS particles
    
    
	//double PTW = 0.;
	double pt_miss = p_tot.polarRadius();
    
    _h_pt_miss->fill (pt_miss, weight);
    /**/ 
    if (pt_miss < 25.) // cut: PT_miss >= 25 GeV
      return;    
    /**/
    
    int nlep = theIsolChargedLeptons.size();
    //if (nlep == 2) { //TODO: what if there are 2 l+-'s not consistent with Z?
    if (nlep == 2 && theIsolChargedLeptons[0].getPdgId() == -theIsolChargedLeptons[1].getPdgId() ) { //two same flavour leptons, opposite charge
      FourMomentum psum;
      for (ParticleVector::const_iterator p = theIsolChargedLeptons.begin(); p != theIsolChargedLeptons.end(); ++p) {
        FourMomentum plep(p->getMomentum());
        psum += plep;
      }
      // Z candidate - drop the event
      if (psum.mass() - 91.18 < 15.0) return; 
	}
	else if (nlep == 1) {
      Particle& p = theIsolChargedLeptons[0];
      FourMomentum plep(p.getMomentum());
      FourMomentum PTW = FourMomentum(plep-p_tot).pz(0.0);
      log << Log::DEBUG << " PTW: " << PTW.polarRadius() << endl;
      _h_pt_e->fill(plep.polarRadius(), weight);
      _h_pt_w->fill(PTW.polarRadius(), weight);
	}
      
    log << Log::DEBUG << "Finished analyzing" << endl;
}


void D0_2001_S4674421::finalize() { 
  //Log& log = getLog();
  _h_w_rec_eff_factor = histogramFactory().divide("/D0_2001_S4674421/W_rec_eff", *_h_pt_w_true, *_h_pt_w); // true/measured(after cuts) => factors large one!
  
  //histogramFactory().destroy(_h_pt_w);
  //histogramFactory().destroy(_h_pt_w_true);
}


