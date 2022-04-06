// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief DELPHI W decay analysis
  class DELPHI_2001_I526164 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_2001_I526164);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // Hisotgram booking
      // qqbar
      // counters for multiplicities
      for(unsigned int ix=0;ix<7;++ix) {
        std::ostringstream title;
        title << "_n_qq_" << ix;
        book(_n_qq[ix],title.str());
      }
      // spectra
      //_WW=false;
      if(isCompatibleWithSqrtS(183.)) {
        book(_h_qq_K0 ,16,1,1);
        book(_h_qq_Lam,16,1,3);
        book(_h_p_charged[0],7,1,2);
        book(_h_p_charged[1],7,1,1);
        book(_h_p_chargedB[0],"/TMP/h_7_1_2",refData(8,1,1));
        book(_h_p_chargedB[1],"/TMP/h_7_1_1",refData(8,1,1));
        book(_h_xi_charged [0],10,1,2);
        book(_h_xi_charged [1],10,1,1);
        book(_h_xi_chargedB[0],"/TMP/h_10_1_2",refData(10,1,3));
        book(_h_xi_chargedB[1],"/TMP/h_10_1_1",refData(10,1,3));
        book(_h_pT_charged [0],12,1,2);
        book(_h_pT_charged [1],12,1,1);
        book(_h_pT_chargedB[0],"/TMP/h_12_1_2",refData(12,1,3));
        book(_h_pT_chargedB[1],"/TMP/h_12_1_1",refData(12,1,3));
        //_WW=true;
      }
      else if (isCompatibleWithSqrtS(189.)) {
        book(_h_qq_K0 ,16,1,2);
        book(_h_qq_Lam,16,1,4);
        book(_h_p_charged [0],5,1,2);
        book(_h_p_charged [1],5,1,1);
        book(_h_p_chargedB[0],"/TMP/h_5_1_2",refData(6,1,1));
        book(_h_p_chargedB[1],"/TMP/h_5_1_1",refData(6,1,1));
        book(_h_xi_charged [0],9,1,2);
        book(_h_xi_charged [1],9,1,1);
        book(_h_xi_chargedB[0],"/TMP/h_9_1_2",refData(9,1,3));
        book(_h_xi_chargedB[1],"/TMP/h_9_1_1",refData(9,1,3));
        book(_h_pT_charged [0],11,1,2);
        book(_h_pT_charged [1],11,1,1);
        book(_h_pT_chargedB[0],"/TMP/h_11_1_2",refData(11,1,3));
        book(_h_pT_chargedB[1],"/TMP/h_11_1_1",refData(11,1,3));
        for(unsigned int ix=0;ix<3;++ix) {
          for(unsigned int iy=0;iy<4;++iy) {
            book(_h_xi_ident[ix][iy],13+ix,1,iy+1);
          }
        }
        //_WW=true;
      }
      /*if(_WW) {
        for(unsigned int ix=0;ix<2;++ix) {
          for(unsigned int iy=0;iy<5;++iy) {
            std::ostringstream title;
            title << "_n_WW_" << ix << "_" << iy << "\n";
            book(_n_WW[ix][iy],title.str());
          }
        }
      }*/
    }

    /*void findChildren(Particles & children, set<long> & barcodes, ConstGenParticlePtr parent) {
      ConstGenVertexPtr dv = parent->end_vertex();
      for(ConstGenParticlePtr child :  HepMCUtils::particles(dv, Relatives::CHILDREN)) {
        long pid = abs(child->pdg_id());
        if(pid==PID::PIPLUS|| pid==PID::PROTON || pid==PID::KPLUS) {
          if(barcodes.find(child->barcode())==barcodes.end()) {
            children.push_back(Particle(child));
            barcodes.insert(child->barcode());
          }
        }
        else if(child->end_vertex()) {
          findChildren(children,barcodes,child);
        }
        else {
          if(PID::isCharged(child->pdg_id())) {
            if(barcodes.find(child->barcode())==barcodes.end()) {
              children.push_back(Particle(child));
              barcodes.insert(child->barcode());
            }
          }
        }
      }
    }*/

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // todo this is not ideal, based on initialquarks projection
      // Classify the event qqbar semi-l W or hadronic W (veto taus)
      // type of W decay true hadronic, false leptonic taus are vetoed
      bool hadronic_w[2]={false,false};
      // type of event
      bool qqbar=false;
      // store the Ws
      ConstGenParticlePtr Wbosons[2];
      // find the type of event
      unsigned int iw=0;
      for (ConstGenParticlePtr p : HepMCUtils::particles(event.genEvent())) {
        ConstGenVertexPtr pv = p->production_vertex();
        ConstGenVertexPtr dv = p->end_vertex();
        const PdgId pid = abs(p->pdg_id());
        bool passed = inRange((long)pid, 1, 6) || pid==24;
	
        if (passed) {
          if (pv != 0) {
            for (ConstGenParticlePtr pp : HepMCUtils::particles(pv, Relatives::PARENTS)){
              // Only accept if parent is electron or Z0 or photon
              const PdgId pid = abs(pp->pdg_id());
              passed = (pid == PID::ELECTRON || ((abs(pp->pdg_id()) == PID::ZBOSON || abs(pp->pdg_id()) == PID::GAMMA)&&pp->generated_mass()>100.));
            }
          }
          else {
            passed = false;
          }
        }
        if(passed) {
          if (dv != 0 && pid==24) {
            ConstGenParticlePtr parent = p;
            std::vector<ConstGenParticlePtr> outgoing = HepMCUtils::particles(dv, Relatives::CHILDREN);
            do {
              if(outgoing.size()==1 && outgoing[0]->pdg_id()==p->pdg_id()) {
                parent = outgoing[0];
                dv = parent->end_vertex();
                outgoing = HepMCUtils::particles(dv, Relatives::CHILDREN);
              }
              else
                break;
            }
            while(true);
            // type of W decay
            Wbosons[iw] = parent;
            for (ConstGenParticlePtr pp : outgoing) {
              // veto taus
              if(abs(pp->pdg_id())==15 || abs(pp->pdg_id())==16)
                vetoEvent;
              else if(abs(pp->pdg_id())==11 || abs(pp->pdg_id())==12 ||
                abs(pp->pdg_id())==13 || abs(pp->pdg_id())==14)
                  hadronic_w[iw]=false;
              else if(abs(pp->pdg_id())<=5)
                hadronic_w[iw]=true;
            }
            ++iw;
          }
          else {
            qqbar=true;
            break;
          }
        }
      }
      // veto leptonic events
      if(!qqbar && !hadronic_w[0] && !hadronic_w[1]) vetoEvent;
      // now do the analysis
      if(qqbar) {
        _n_qq[0]->fill();
        UnstableParticles ufs=apply<UnstableParticles>(event,"UFS");
        for(const Particle & p : ufs.particles()) {
          if(p.abspid()==PID::PIPLUS) {
            _n_qq[1]->fill();
            _n_qq[2]->fill();
          }
          else if(p.abspid()==PID::KPLUS) {
            _n_qq[1]->fill();
            _n_qq[3]->fill();
          }
          else if(p.abspid()==PID::PROTON) {
            _n_qq[1]->fill();
            _n_qq[5]->fill();
          }
          else if(p.abspid()==130 || p.abspid()==310 ) {
            _n_qq[4]->fill();
            if(_h_qq_K0) _h_qq_K0->fill(-log(2.*p.p3().mod()/sqrtS()));
          }
          else if(p.abspid()==PID::LAMBDA ) {
            _n_qq[6]->fill();
            if(_h_qq_Lam) _h_qq_Lam->fill(-log(2.*p.p3().mod()/sqrtS()));
          }
        }
      }
      // WW production
      /*else if(_WW) {
        const double mW=80.379*GeV;
        // find W decay products
        Particles particles;
        set<long> barcodes;
        for(unsigned int ix=0;ix<2;++ix)
          if(hadronic_w[ix]) findChildren(particles,barcodes,Wbosons[ix]);
        // calculate thrust
        Thrust thrust;
        thrust.calc(particles);
        // type of event
        unsigned int itype= !hadronic_w[0]||!hadronic_w[1] ? 0 : 1;
        _n_WW[itype][0]->fill();
        // loop over particles and fill histos
        for(const Particle & p : particles) {
          // Get momentum of each particle.
          const Vector3 mom3 = p.p3();
          // Scaled momenta.
          const double mom = mom3.mod();
          const double xiW = -log(2.*mom/mW);
          const double xiL = -log(2.*mom/sqrtS());
          // Get momenta components w.r.t. thrust
          const double pTinT = dot(mom3, thrust.thrustMajorAxis());
          const double pToutT = dot(mom3, thrust.thrustMinorAxis());
          double pT = sqrt(sqr(pTinT)+sqr(pToutT));
          // fill charged particle hists
          _h_p_charged  [itype]->fill(mom);
          _h_xi_charged [itype]->fill(xiL);
          _h_pT_charged [itype]->fill(pT );
          _h_p_chargedB [itype]->fill(mom);
          _h_xi_chargedB[itype]->fill(xiL);
          _h_pT_chargedB[itype]->fill(pT );
          // and identified particles
          if(_h_xi_ident[itype][0]) _h_xi_ident[itype][0]->fill(xiW);
          _n_WW[itype][1]->fill();
          if(p.abspid()==PID::PIPLUS) {
            if(_h_xi_ident[itype][1]) _h_xi_ident[itype][1]->fill(xiW);
            _n_WW[itype][2]->fill();
          }
          else if(p.abspid()==PID::KPLUS) {
            if(_h_xi_ident[itype][2]) _h_xi_ident[itype][2]->fill(xiW);
            _n_WW[itype][3]->fill();
          }
          else if(p.abspid()==PID::PROTON){
            if(_h_xi_ident[itype][3]) _h_xi_ident[itype][3]->fill(xiW);
            _n_WW[itype][4]->fill();
          }
        }
        // boosted specta in W rest frame
        if(itype==0) {
          // boost to rest frame
          Particle Whad = Particle( hadronic_w[0] ? Wbosons[0] : Wbosons[1]); 
          LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(Whad.momentum().betaVec());
          // spectra
          for(const Particle & p : particles) {
            // Get momentum of each particle.
            FourMomentum phad = boost.transform(p.momentum());
            const Vector3 mom3 = phad.p3();
            // Scaled momenta.
            const double mom = mom3.mod();
            const double scaledMom = 2.*mom/mW;
            const double xi = -log(scaledMom);
            // /identified particle spectra
            if(_h_xi_ident[2][0]) _h_xi_ident[2][0]->fill(xi);
            if(p.abspid()==PID::PIPLUS) {
              if(_h_xi_ident[2][1]) _h_xi_ident[2][1]->fill(xi);
            }
            else if(p.abspid()==PID::KPLUS) {
              if(_h_xi_ident[2][2]) _h_xi_ident[2][2]->fill(xi);
            }
            else if(p.abspid()==PID::PROTON){
              if(_h_xi_ident[2][3]) _h_xi_ident[2][3]->fill(xi);
            }
          }
        }
      }*/
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // qq
      if(_n_qq[0]->effNumEntries()!=0.) {
        // scale spectra if needed
        if(_h_qq_K0 ) scale(_h_qq_K0 , 1./ *_n_qq[0]);
        if(_h_qq_Lam) scale(_h_qq_Lam, 1./ *_n_qq[0]);
        for(unsigned int ih=1;ih<7;++ih) {
          if(_n_qq[ih]->effNumEntries()==0.) continue;
          unsigned int ix= ih==1 ? 1 : 2, iy= ih==1 ? 1 : ih-1;
          // ratio
          std::ostringstream title;
          title << "/TMP/n_" << ix << "_" << iy;
          Scatter1D sTemp(title.str());
          Scatter1DPtr s1d = registerAO(sTemp);
          divide(_n_qq[ih],_n_qq[0],s1d);
          // hist for axis
          Scatter2D temphisto(refData(ix, 1, iy));
          Scatter2DPtr ratio;
          book(ratio, ix, 1, iy);
          for (size_t b = 0; b < temphisto.numPoints(); b++) {
            const double x  = temphisto.point(b).x();
            pair<double,double> ex = temphisto.point(b).xErrs();
            pair<double,double> ex2 = ex;
            if(ex2.first ==0.) ex2. first=0.0001;
            if(ex2.second==0.) ex2.second=0.0001;
            if (inRange(sqrtS(), x-ex2.first, x+ex2.second)) {
              ratio->addPoint(x,s1d->points()[0].x(),ex,s1d->points()[0].xErrs());
            }
            else {
              ratio->addPoint(x, 0., ex, make_pair(0.,.0));
            }
          }
        }
      }
      // WW
      /*if( _WW && 
	  (_n_WW[0][0]->effNumEntries()!=0. ||
	   _n_WW[1][0]->effNumEntries()!=0. )) {
        Scatter1D ratios[2];
        unsigned int iloc = isCompatibleWithSqrtS(183.) ? 1 : 0;
        for(unsigned int iw=0;iw<2;++iw) {
          if(_n_WW[iw][0]->effNumEntries()==0.) continue;
          // scale distributions
          scale(_h_p_charged[iw] , 1./ *_n_WW[iw][0]);
          scale(_h_xi_charged[iw], 1./ *_n_WW[iw][0]);
          scale(_h_pT_charged[iw], 1./ *_n_WW[iw][0]);
          scale(_h_p_chargedB[iw] , 1./ *_n_WW[iw][0]);
          scale(_h_xi_chargedB[iw], 1./ *_n_WW[iw][0]);
          scale(_h_pT_chargedB[iw], 1./ *_n_WW[iw][0]);
          for(unsigned int iy=0;iy<4;++iy) {
            if(_h_xi_ident[iw][iy]) {
              scale(_h_xi_ident[iw][iy], 1./ *_n_WW[iw][0]);
              if(iw==0)
                scale(_h_xi_ident[2][iy], 1./ *_n_WW[iw][0]);
            }
          }
          // charge mult
          ratios[iw] =  *_n_WW[iw][1]/ *_n_WW[iw][0];
          // sort out identified particle multiplicities
          if(iloc==0) {
            for(unsigned int iy=2;iy<5;++iy) {
              Scatter1D R = *_n_WW[iw][iy]/ *_n_WW[iw][0];
              Scatter2D temphisto(refData(4, 1, 2*(iy-1)-iw));
              Scatter2DPtr mult;
              book(mult, 4, 1,  2*(iy-1)-iw);
              for (size_t b = 0; b < temphisto.numPoints(); b++) {
                const double x  = temphisto.point(b).x();
                pair<double,double> ex = temphisto.point(b).xErrs();
                pair<double,double> ex2 = ex;
                if(ex2.first ==0.) ex2. first=0.0001;
                if(ex2.second==0.) ex2.second=0.0001;
                if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
                  mult->addPoint(x, R.point(0).x(), ex, R.point(0).xErrs());
                }
                else {
                  mult->addPoint(x, 0., ex, make_pair(0.,.0));
                }
              }
            }
          }
        }
        // charged mults
        Scatter2D temphisto(refData(3, 1, 1));
        Scatter2DPtr mult;
        book(mult,3,1,1);
        for (size_t b = 0; b < temphisto.numPoints(); b++) {
          const double x  = temphisto.point(b).x();
          pair<double,double> ex = temphisto.point(b).xErrs();
          if((iloc==1 && b==0) || (iloc==0 && b==1))
            mult->addPoint(x, ratios[1].point(0).x(), ex, ratios[1].point(0).xErrs());
          else if((iloc==1 && b==2) || (iloc==0 && b==3))
            mult->addPoint(x, ratios[0].point(0).x(), ex, ratios[0].point(0).xErrs());
          else
            mult->addPoint(x, 0., ex, make_pair(0.,.0));
        }
        // difference histos
        // momentum
        Scatter2D htemp = mkScatter(*(_h_p_chargedB[0]));
        htemp.scaleY(2.);
        Scatter2DPtr hnew;
        book(hnew,6+2*iloc,1,1);
        *hnew = YODA::subtract(*(_h_p_chargedB[1]),htemp);
        hnew->setPath("/"+name()+"/"+mkAxisCode(6+2*iloc,1,1));
        // xi
        Scatter2D htemp2 = mkScatter(*(_h_xi_chargedB[0]));
        htemp2.scaleY(2.);
        Scatter2DPtr hnew2;
        book(hnew2,9+iloc,1,3);
        *hnew2 = YODA::subtract(*(_h_xi_chargedB[1]),htemp2);
        hnew2->setPath("/"+name()+"/"+mkAxisCode(9+iloc,1,3));
        // pt
        Scatter2D hytemp3 = mkScatter(*(_h_pT_chargedB[0]));
        hytemp3.scaleY(2.);
        Scatter2DPtr hnew3;
        book(hnew3,11+iloc,1,3);
        *hnew3 = YODA::subtract(*(_h_pT_chargedB[1]),hytemp3);
        hnew3->setPath("/"+name()+"/"+mkAxisCode(11+iloc,1,3));
      }*/
    }

    ///@}


    /// @name Histograms
    ///@{
    // qqbar
    Histo1DPtr _h_qq_K0,_h_qq_Lam;
    CounterPtr _n_qq[7];
    // WW
    Histo1DPtr _h_p_charged [2],_h_xi_charged [2],_h_pT_charged [2],_h_xi_ident[3][4];
    Histo1DPtr _h_p_chargedB[2],_h_xi_chargedB[2],_h_pT_chargedB[2];
    //CounterPtr _n_WW[2][5];
    //bool _WW;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(DELPHI_2001_I526164);

}
