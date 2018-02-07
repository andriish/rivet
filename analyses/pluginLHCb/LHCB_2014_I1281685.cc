// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Charged particle multiplicities and densities in $pp$ collisions at $\sqrt{s} = 7$ TeV
  class LHCB_2014_I1281685 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    LHCB_2014_I1281685()
      : Analysis("LHCB_2014_I1281685"),
        _p_min(2.0),
        _pt_min(0.2),
        _eta_min(2.0),
        _eta_max(4.8),
        _maxlft(1.0e-11)
    {    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      fillMap(_partLftMap);

      // Projections
      declare(ChargedFinalState(_eta_min, _eta_max, _pt_min*GeV), "CFS");

      // Book histograms
      book(_h_mult_total  ,"d03-x01-y01", 50, 0.5, 50.5);

      book(_h_mult_eta[0] ,"d04-x01-y01", 21, -0.5, 20.5); //eta=[2.0,2.5]
      book(_h_mult_eta[1] ,"d04-x01-y02", 21, -0.5, 20.5); //eta=[2.5,3.0]
      book(_h_mult_eta[2] ,"d04-x01-y03", 21, -0.5, 20.5); //eta=[3.0,3.5]
      book(_h_mult_eta[3] ,"d04-x01-y04", 21, -0.5, 20.5); //eta=[3.5,4.0]
      book(_h_mult_eta[4] ,"d04-x01-y05", 21, -0.5, 20.5); //eta=[4.0,4.5]

      book(_h_mult_pt[0]  ,"d05-x01-y01", 21, -0.5, 20.5); //pT=[0.2,0.3]GeV
      book(_h_mult_pt[1]  ,"d05-x01-y02", 21, -0.5, 20.5); //pT=[0.3,0.4]GeV
      book(_h_mult_pt[2]  ,"d05-x01-y03", 21, -0.5, 20.5); //pT=[0.4,0.6]GeV
      book(_h_mult_pt[3]  ,"d05-x01-y04", 21, -0.5, 20.5); //pT=[0.6,1.0]GeV
      book(_h_mult_pt[4]  ,"d05-x01-y05", 21, -0.5, 20.5); //pT=[1.0,2.0]GeV

      book(_h_dndeta      ,"d01-x01-y01", 14, 2.0, 4.8); //eta=[2,4.8]
      book(_h_dndpt       ,"d02-x01-y01", 18, 0.2, 2.0); //pT =[0,2]GeV

      // Counters
      book(_sumW, "TMP/sumW");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Variable to store multiplicities per event
      int LHCbcountAll = 0; //count particles fulfiling all requirements
      int LHCbcountEta[8] = {0,0,0,0,0,0,0,0}; //count per eta-bin
      int LHCbcountPt[7]  = {0,0,0,0,0,0,0};   //count per pT-bin
      vector<double> val_dNdEta;
      vector<double> val_dNdPt;
      val_dNdEta.clear();
      val_dNdPt.clear();

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      foreach (const Particle& p, cfs.particles()) {
        int id = p.pdgId();
        // continue if particle is not a pion, kaon, proton, muon or electron
        if ( !( (abs(id) == 211) || (abs(id) == 321) || (abs(id) == 2212) || (abs(id) == 13) || (abs(id) == 11)) ) {
          continue;
        }

        const FourMomentum& qmom = p.momentum();
        const double eta = p.momentum().eta();
        const double pT  = p.momentum().pT();
        //minimum momentum
        if (qmom.p3().mod() < _p_min) continue;
        //minimum tr. momentum
        if (pT < _pt_min) continue;
        //eta range
        if ((eta < _eta_min) || (eta > _eta_max)) continue;

        /* Select only prompt particles via lifetime */
        //Sum of all mother lifetimes (PDG lifetime) < 10ps
        double ancestors_sumlft = getAncestorSumLifetime(p);
        if( (ancestors_sumlft > _maxlft) || (ancestors_sumlft < 0) ) continue;

        //after all cuts;
        LHCbcountAll++; //count particles in whole kin. range

        //in eta bins
        if( eta >2.0 && eta <= 2.5) LHCbcountEta[0]++;
        if( eta >2.5 && eta <= 3.0) LHCbcountEta[1]++;
        if( eta >3.0 && eta <= 3.5) LHCbcountEta[2]++;
        if( eta >3.5 && eta <= 4.0) LHCbcountEta[3]++;
        if( eta >4.0 && eta <= 4.5) LHCbcountEta[4]++;
        if( eta >2.0 && eta <= 4.8) LHCbcountEta[5]++; //cross-check
        //in pT bins
        if( pT > 0.2 && pT <= 0.3) LHCbcountPt[0]++;
        if( pT > 0.3 && pT <= 0.4) LHCbcountPt[1]++;
        if( pT > 0.4 && pT <= 0.6) LHCbcountPt[2]++;
        if( pT > 0.6 && pT <= 1.0) LHCbcountPt[3]++;
        if( pT > 1.0 && pT <= 2.0) LHCbcountPt[4]++;
        if( pT > 0.2)              LHCbcountPt[5]++;   //cross-check

        //particle densities -> need proper normalization (finalize)
        val_dNdPt.push_back( pT );
        val_dNdEta.push_back( eta );
      }//end foreach


      // Fill histograms only, if at least 1 particle pre event was within the
      // kinematic range of the analysis!
      if (LHCbcountAll) {
        _sumW->fill();

        _h_mult_total->fill(LHCbcountAll);

        _h_mult_eta[0]->fill(LHCbcountEta[0]);
        _h_mult_eta[1]->fill(LHCbcountEta[1]);
        _h_mult_eta[2]->fill(LHCbcountEta[2]);
        _h_mult_eta[3]->fill(LHCbcountEta[3]);
        _h_mult_eta[4]->fill(LHCbcountEta[4]);

        _h_mult_pt[0]->fill(LHCbcountPt[0]);
        _h_mult_pt[1]->fill(LHCbcountPt[1]);
        _h_mult_pt[2]->fill(LHCbcountPt[2]);
        _h_mult_pt[3]->fill(LHCbcountPt[3]);
        _h_mult_pt[4]->fill(LHCbcountPt[4]);

        for (size_t part = 0; part < val_dNdEta.size(); part++)
          _h_dndeta->fill(val_dNdEta[part]);

        for (size_t part = 0; part < val_dNdPt.size(); part++)
          _h_dndpt->fill(val_dNdPt[part]);

      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double scalefactor = 1.0/_sumW->val(); // normalize multiplicity histograms by nEvents
      const double scale1k = 1000.; // to match '10^3' scale in reference histograms

      scale( _h_dndeta, scalefactor );
      scale( _h_dndpt,  scalefactor*0.1  ); //additional factor 0.1 for [0.1 GeV/c]
      scale( _h_mult_total, scalefactor*scale1k);

      _h_mult_eta[0]->scaleW( scalefactor*scale1k );
      _h_mult_eta[1]->scaleW( scalefactor*scale1k );
      _h_mult_eta[2]->scaleW( scalefactor*scale1k );
      _h_mult_eta[3]->scaleW( scalefactor*scale1k );
      _h_mult_eta[4]->scaleW( scalefactor*scale1k );

      _h_mult_pt[0]->scaleW( scalefactor*scale1k );
      _h_mult_pt[1]->scaleW( scalefactor*scale1k );
      _h_mult_pt[2]->scaleW( scalefactor*scale1k );
      _h_mult_pt[3]->scaleW( scalefactor*scale1k );
      _h_mult_pt[4]->scaleW( scalefactor*scale1k );
    }

    //@}


  private:


    // Get mean PDG lifetime for particle with PID
    double getLifetime(int pid) {
      double lft = 0.;
      map<int, double>::iterator pPartLft = _partLftMap.find(pid);
      if (pPartLft != _partLftMap.end()) {
        lft = (*pPartLft).second;
      } else {
        // allow identifying missing life times only in debug mode
        MSG_DEBUG("Could not determine lifetime for particle with PID " << pid  << "... Assume non-prompt particle");
        lft = -1;
      }
      return lft;
    }


    // Get sum of all ancestor particles
    const double getAncestorSumLifetime(const Particle& p) {
      double lftSum = 0.;
      double plft = 0.;
      const GenParticle* part = p.genParticle();
      if ( 0 == part ) return -1;
      const GenVertex* ivtx = part->production_vertex();
      while(ivtx) {
        if (ivtx->particles_in_size() < 1) { lftSum = -1.; break; };
        const GenVertex::particles_in_const_iterator iPart_invtx = ivtx->particles_in_const_begin();
        part = (*iPart_invtx);
        if ( !(part) ) { lftSum = -1.; break; };
        ivtx = part->production_vertex();
        if ( (part->pdg_id() == 2212) || !(ivtx) ) break; // reached beam
        plft = getLifetime(part->pdg_id());
        if (plft < 0.) { lftSum = -1.; break; };
        lftSum += plft;
      }
      return (lftSum);
    }


    /// Hard-coded map linking PDG ID with PDG lifetime[s] (converted from ParticleTable.txt)
    bool fillMap(map<int, double>& m) {
      // PDGID  = LIFETIME
      m[22] = 1.000000e+016;
      m[-11] = 1.000000e+016;
      m[11] = 1.000000e+016;
      m[12] = 1.000000e+016;
      m[-13] = 2.197036e-006;
      m[13] = 2.197036e-006;
      m[111] = 8.438618e-017;
      m[211] = 2.603276e-008;
      m[-211] = 2.603276e-008;
      m[130] = 5.174624e-008;
      m[321] = 1.238405e-008;
      m[-321] = 1.238405e-008;
      m[2112] =    885.646128;
      m[2212] = 1.000000e+016;
      m[-2212] = 1.000000e+016;
      m[310] = 8.934603e-011;
      m[221] = 5.578070e-019;
      m[3122] = 2.631796e-010;
      m[3222] = 8.018178e-011;
      m[3212] = 7.395643e-020;
      m[3112] = 1.479129e-010;
      m[3322] = 2.899613e-010;
      m[3312] = 1.637344e-010;
      m[3334] = 8.207135e-011;
      m[-2112] =    885.646128;
      m[-3122] = 2.631796e-010;
      m[-3222] = 8.018178e-011;
      m[-3212] = 7.395643e-020;
      m[-3112] = 1.479129e-010;
      m[-3322] = 2.899613e-010;
      m[-3312] = 1.637344e-010;
      m[-3334] = 8.207135e-011;
      m[113] = 4.411610e-024;
      m[213] = 4.411610e-024;
      m[-213] = 4.411610e-024;
      m[223] = 7.798723e-023;
      m[333] = 1.545099e-022;
      m[323] = 1.295693e-023;
      m[-323] = 1.295693e-023;
      m[313] = 1.298249e-023;
      m[-313] = 1.298249e-023;
      m[20213] = 1.500000e-024;
      m[-20213] = 1.500000e-024;
      m[450000000] = 1.000000e+015;
      m[460000000] = 1.000000e+015;
      m[470000000] = 1.000000e+015;
      m[480000000] = 1.000000e+015;
      m[490000000] = 1.000000e+015;
      m[20022] = 1.000000e+016;
      m[-15] = 2.906014e-013;
      m[15] = 2.906014e-013;
      m[24] = 3.104775e-025;
      m[-24] = 3.104775e-025;
      m[23] = 2.637914e-025;
      m[411] = 1.051457e-012;
      m[-411] = 1.051457e-012;
      m[421] = 4.116399e-013;
      m[-421] = 4.116399e-013;
      m[431] = 4.904711e-013;
      m[-431] = 4.904711e-013;
      m[4122] = 1.994582e-013;
      m[-4122] = 1.994582e-013;
      m[443] = 7.565657e-021;
      m[413] = 6.856377e-021;
      m[-413] = 6.856377e-021;
      m[423] = 1.000003e-019;
      m[-423] = 1.000003e-019;
      m[433] = 1.000003e-019;
      m[-433] = 1.000003e-019;
      m[521] = 1.671000e-012;
      m[-521] = 1.671000e-012;
      m[511] = 1.536000e-012;
      m[-511] = 1.536000e-012;
      m[531] = 1.461000e-012;
      m[-531] = 1.461000e-012;
      m[541] = 4.600000e-013;
      m[-541] = 4.600000e-013;
      m[5122] = 1.229000e-012;
      m[-5122] = 1.229000e-012;
      m[4112] = 4.388081e-022;
      m[-4112] = 4.388081e-022;
      m[4212] = 3.999999e-022;
      m[-4212] = 3.999999e-022;
      m[4222] = 3.291060e-022;
      m[-4222] = 3.291060e-022;
      m[25] = 9.400000e-026;
      m[35] = 9.400000e-026;
      m[36] = 9.400000e-026;
      m[37] = 9.400000e-026;
      m[-37] = 9.400000e-026;
      m[4312] = 9.800002e-014;
      m[-4312] = 9.800002e-014;
      m[4322] = 3.500001e-013;
      m[-4322] = 3.500001e-013;
      m[4332] = 6.453061e-014;
      m[-4332] = 6.453061e-014;
      m[4132] = 9.824063e-014;
      m[-4132] = 9.824063e-014;
      m[4232] = 4.417532e-013;
      m[-4232] = 4.417532e-013;
      m[5222] = 1.000000e-019;
      m[-5222] = 1.000000e-019;
      m[5212] = 1.000000e-019;
      m[-5212] = 1.000000e-019;
      m[5112] = 1.000000e-019;
      m[-5112] = 1.000000e-019;
      m[5312] = 1.000000e-019;
      m[-5312] = 1.000000e-019;
      m[5322] = 1.000000e-019;
      m[-5322] = 1.000000e-019;
      m[5332] = 1.550000e-012;
      m[-5332] = 1.550000e-012;
      m[5132] = 1.390000e-012;
      m[-5132] = 1.390000e-012;
      m[5232] = 1.390000e-012;
      m[-5232] = 1.390000e-012;
      m[100443] = 2.194041e-021;
      m[331] = 3.258476e-021;
      m[441] = 4.113826e-023;
      m[10441] = 4.063038e-023;
      m[20443] = 7.154480e-022;
      m[445] = 3.164482e-022;
      m[9000111] = 1.149997e-023;
      m[9000211] = 1.149997e-023;
      m[-9000211] = 1.149997e-023;
      m[20113] = 1.500000e-024;
      m[115] = 6.151516e-024;
      m[215] = 6.151516e-024;
      m[-215] = 6.151516e-024;
      m[10323] = 7.313469e-024;
      m[-10323] = 7.313469e-024;
      m[10313] = 7.313469e-024;
      m[-10313] = 7.313469e-024;
      m[20323] = 3.782829e-024;
      m[-20323] = 3.782829e-024;
      m[20313] = 3.782829e-024;
      m[-20313] = 3.782829e-024;
      m[10321] = 2.238817e-024;
      m[-10321] = 2.238817e-024;
      m[10311] = 2.238817e-024;
      m[-10311] = 2.238817e-024;
      m[325] = 6.682357e-024;
      m[-325] = 6.682357e-024;
      m[315] = 6.038644e-024;
      m[-315] = 6.038644e-024;
      m[10411] = 4.380000e-024;
      m[20413] = 2.630000e-024;
      m[10413] = 3.290000e-023;
      m[-415] = 2.632849e-023;
      m[-10411] = 4.380000e-024;
      m[-20413] = 2.630000e-024;
      m[-10413] = 3.290000e-023;
      m[415] = 2.632849e-023;
      m[10421] = 4.380000e-024;
      m[20423] = 2.630000e-024;
      m[10423] = 3.482604e-023;
      m[-425] = 2.861792e-023;
      m[-10421] = 4.380000e-024;
      m[-20423] = 2.630000e-024;
      m[-10423] = 3.482604e-023;
      m[425] = 2.861792e-023;
      m[10431] = 6.582100e-022;
      m[20433] = 6.582100e-022;
      m[10433] = 6.582100e-022;
      m[435] = 4.388100e-023;
      m[-10431] = 6.582100e-022;
      m[-20433] = 6.582100e-022;
      m[-10433] = 6.582100e-022;
      m[-435] = 4.388100e-023;
      m[2224] = 5.485102e-024;
      m[2214] = 5.485102e-024;
      m[2114] = 5.485102e-024;
      m[1114] = 5.485102e-024;
      m[-2224] = 5.485102e-024;
      m[-2214] = 5.485102e-024;
      m[-2114] = 5.485102e-024;
      m[-1114] = 5.485102e-024;
      m[-523] = 1.000019e-019;
      m[523] = 1.000019e-019;
      m[513] = 1.000019e-019;
      m[-513] = 1.000019e-019;
      m[533] = 1.000000e-019;
      m[-533] = 1.000000e-019;
      m[10521] = 4.390000e-024;
      m[20523] = 2.630000e-024;
      m[10523] = 1.650000e-023;
      m[525] = 1.310000e-023;
      m[-10521] = 4.390000e-024;
      m[-20523] = 2.630000e-024;
      m[-10523] = 1.650000e-023;
      m[-525] = 1.310000e-023;
      m[10511] = 4.390000e-024;
      m[20513] = 2.630000e-024;
      m[10513] = 1.650000e-023;
      m[515] = 1.310000e-023;
      m[-10511] = 4.390000e-024;
      m[-20513] = 2.630000e-024;
      m[-10513] = 1.650000e-023;
      m[-515] = 1.310000e-023;
      m[10531] = 4.390000e-024;
      m[20533] = 2.630000e-024;
      m[10533] = 1.650000e-023;
      m[535] = 1.310000e-023;
      m[-10531] = 4.390000e-024;
      m[-20533] = 2.630000e-024;
      m[-10533] = 1.650000e-023;
      m[-535] = 1.310000e-023;
      m[14] = 1.000000e+016;
      m[-14] = 1.000000e+016;
      m[-12] = 1.000000e+016;
      m[1] = 0.000000e+000;
      m[-1] = 0.000000e+000;
      m[2] = 0.000000e+000;
      m[-2] = 0.000000e+000;
      m[3] = 0.000000e+000;
      m[-3] = 0.000000e+000;
      m[4] = 0.000000e+000;
      m[-4] = 0.000000e+000;
      m[5] = 0.000000e+000;
      m[-5] = 0.000000e+000;
      m[6] = 4.707703e-025;
      m[-6] = 4.707703e-025;
      m[7] = 0.000000e+000;
      m[-7] = 0.000000e+000;
      m[8] = 0.000000e+000;
      m[-8] = 0.000000e+000;
      m[16] = 1.000000e+016;
      m[-16] = 1.000000e+016;
      m[17] = 0.000000e+000;
      m[-17] = 0.000000e+000;
      m[18] = 0.000000e+000;
      m[-18] = 0.000000e+000;
      m[21] = 0.000000e+000;
      m[32] = 0.000000e+000;
      m[33] = 0.000000e+000;
      m[34] = 0.000000e+000;
      m[-34] = 0.000000e+000;
      m[39] = 0.000000e+000;
      m[41] = 0.000000e+000;
      m[-41] = 0.000000e+000;
      m[42] = 0.000000e+000;
      m[-42] = 0.000000e+000;
      m[43] = 0.000000e+000;
      m[44] = 0.000000e+000;
      m[-44] = 0.000000e+000;
      m[81] = 0.000000e+000;
      m[82] = 0.000000e+000;
      m[-82] = 0.000000e+000;
      m[83] = 0.000000e+000;
      m[84] = 3.335641e-013;
      m[-84] = 3.335641e-013;
      m[85] = 1.290893e-012;
      m[-85] = 1.290893e-012;
      m[86] = 0.000000e+000;
      m[-86] = 0.000000e+000;
      m[87] = 0.000000e+000;
      m[-87] = 0.000000e+000;
      m[88] = 0.000000e+000;
      m[90] = 0.000000e+000;
      m[91] = 0.000000e+000;
      m[92] = 0.000000e+000;
      m[93] = 0.000000e+000;
      m[94] = 0.000000e+000;
      m[95] = 0.000000e+000;
      m[96] = 0.000000e+000;
      m[97] = 0.000000e+000;
      m[98] = 0.000000e+000;
      m[99] = 0.000000e+000;
      m[117] = 4.088275e-024;
      m[119] = 1.828367e-024;
      m[217] = 4.088275e-024;
      m[-217] = 4.088275e-024;
      m[219] = 1.828367e-024;
      m[-219] = 1.828367e-024;
      m[225] = 3.555982e-024;
      m[227] = 3.917930e-024;
      m[229] = 3.392846e-024;
      m[311] = 1.000000e+016;
      m[-311] = 1.000000e+016;
      m[317] = 4.139699e-024;
      m[-317] = 4.139699e-024;
      m[319] = 3.324304e-024;
      m[-319] = 3.324304e-024;
      m[327] = 4.139699e-024;
      m[-327] = 4.139699e-024;
      m[329] = 3.324304e-024;
      m[-329] = 3.324304e-024;
      m[335] = 8.660687e-024;
      m[337] = 7.565657e-024;
      m[543] = 0.000000e+000;
      m[-543] = 0.000000e+000;
      m[545] = 0.000000e+000;
      m[-545] = 0.000000e+000;
      m[551] = 0.000000e+000;
      m[553] = 1.253738e-020;
      m[555] = 1.000000e+016;
      m[557] = 0.000000e+000;
      m[-450000000] = 0.000000e+000;
      m[-490000000] = 0.000000e+000;
      m[-460000000] = 0.000000e+000;
      m[-470000000] = 0.000000e+000;
      m[1103] = 0.000000e+000;
      m[-1103] = 0.000000e+000;
      m[1112] = 4.388081e-024;
      m[-1112] = 4.388081e-024;
      m[1116] = 1.880606e-024;
      m[-1116] = 1.880606e-024;
      m[1118] = 2.194041e-024;
      m[-1118] = 2.194041e-024;
      m[1212] = 4.388081e-024;
      m[-1212] = 4.388081e-024;
      m[1214] = 5.485102e-024;
      m[-1214] = 5.485102e-024;
      m[1216] = 1.880606e-024;
      m[-1216] = 1.880606e-024;
      m[1218] = 1.462694e-024;
      m[-1218] = 1.462694e-024;
      m[2101] = 0.000000e+000;
      m[-2101] = 0.000000e+000;
      m[2103] = 0.000000e+000;
      m[-2103] = 0.000000e+000;
      m[2116] = 4.388081e-024;
      m[-2116] = 4.388081e-024;
      m[2118] = 2.194041e-024;
      m[-2118] = 2.194041e-024;
      m[2122] = 4.388081e-024;
      m[-2122] = 4.388081e-024;
      m[2124] = 5.485102e-024;
      m[-2124] = 5.485102e-024;
      m[2126] = 1.880606e-024;
      m[-2126] = 1.880606e-024;
      m[2128] = 1.462694e-024;
      m[-2128] = 1.462694e-024;
      m[2203] = 0.000000e+000;
      m[-2203] = 0.000000e+000;
      m[2216] = 4.388081e-024;
      m[-2216] = 4.388081e-024;
      m[2218] = 2.194041e-024;
      m[-2218] = 2.194041e-024;
      m[2222] = 4.388081e-024;
      m[-2222] = 4.388081e-024;
      m[2226] = 1.880606e-024;
      m[-2226] = 1.880606e-024;
      m[2228] = 2.194041e-024;
      m[-2228] = 2.194041e-024;
      m[3101] = 0.000000e+000;
      m[-3101] = 0.000000e+000;
      m[3103] = 0.000000e+000;
      m[-3103] = 0.000000e+000;
      m[3114] = 1.670589e-023;
      m[-3114] = 1.670589e-023;
      m[3116] = 5.485102e-024;
      m[-3116] = 5.485102e-024;
      m[3118] = 3.656734e-024;
      m[-3118] = 3.656734e-024;
      m[3124] = 4.219309e-023;
      m[-3124] = 4.219309e-023;
      m[3126] = 8.227653e-024;
      m[-3126] = 8.227653e-024;
      m[3128] = 3.291061e-024;
      m[-3128] = 3.291061e-024;
      m[3201] = 0.000000e+000;
      m[-3201] = 0.000000e+000;
      m[3203] = 0.000000e+000;
      m[-3203] = 0.000000e+000;
      m[3214] = 1.828367e-023;
      m[-3214] = 1.828367e-023;
      m[3216] = 5.485102e-024;
      m[-3216] = 5.485102e-024;
      m[3218] = 3.656734e-024;
      m[-3218] = 3.656734e-024;
      m[3224] = 1.838582e-023;
      m[-3224] = 1.838582e-023;
      m[3226] = 5.485102e-024;
      m[-3226] = 5.485102e-024;
      m[3228] = 3.656734e-024;
      m[-3228] = 3.656734e-024;
      m[3303] = 0.000000e+000;
      m[-3303] = 0.000000e+000;
      m[3314] = 6.648608e-023;
      m[-3314] = 6.648608e-023;
      m[3324] = 7.233101e-023;
      m[-3324] = 7.233101e-023;
      m[4101] = 0.000000e+000;
      m[-4101] = 0.000000e+000;
      m[4103] = 0.000000e+000;
      m[-4103] = 0.000000e+000;
      m[4114] = 0.000000e+000;
      m[-4114] = 0.000000e+000;
      m[4201] = 0.000000e+000;
      m[-4201] = 0.000000e+000;
      m[4203] = 0.000000e+000;
      m[-4203] = 0.000000e+000;
      m[4214] = 3.291061e-022;
      m[-4214] = 3.291061e-022;
      m[4224] = 0.000000e+000;
      m[-4224] = 0.000000e+000;
      m[4301] = 0.000000e+000;
      m[-4301] = 0.000000e+000;
      m[4303] = 0.000000e+000;
      m[-4303] = 0.000000e+000;
      m[4314] = 0.000000e+000;
      m[-4314] = 0.000000e+000;
      m[4324] = 0.000000e+000;
      m[-4324] = 0.000000e+000;
      m[4334] = 0.000000e+000;
      m[-4334] = 0.000000e+000;
      m[4403] = 0.000000e+000;
      m[-4403] = 0.000000e+000;
      m[4412] = 3.335641e-013;
      m[-4412] = 3.335641e-013;
      m[4414] = 3.335641e-013;
      m[-4414] = 3.335641e-013;
      m[4422] = 3.335641e-013;
      m[-4422] = 3.335641e-013;
      m[4424] = 3.335641e-013;
      m[-4424] = 3.335641e-013;
      m[4432] = 3.335641e-013;
      m[-4432] = 3.335641e-013;
      m[4434] = 3.335641e-013;
      m[-4434] = 3.335641e-013;
      m[4444] = 3.335641e-013;
      m[-4444] = 3.335641e-013;
      m[5101] = 0.000000e+000;
      m[-5101] = 0.000000e+000;
      m[5103] = 0.000000e+000;
      m[-5103] = 0.000000e+000;
      m[5114] = 0.000000e+000;
      m[-5114] = 0.000000e+000;
      m[5142] = 1.290893e-012;
      m[-5142] = 1.290893e-012;
      m[5201] = 0.000000e+000;
      m[-5201] = 0.000000e+000;
      m[5203] = 0.000000e+000;
      m[-5203] = 0.000000e+000;
      m[5214] = 0.000000e+000;
      m[-5214] = 0.000000e+000;
      m[5224] = 0.000000e+000;
      m[-5224] = 0.000000e+000;
      m[5242] = 1.290893e-012;
      m[-5242] = 1.290893e-012;
      m[5301] = 0.000000e+000;
      m[-5301] = 0.000000e+000;
      m[5303] = 0.000000e+000;
      m[-5303] = 0.000000e+000;
      m[5314] = 0.000000e+000;
      m[-5314] = 0.000000e+000;
      m[5324] = 0.000000e+000;
      m[-5324] = 0.000000e+000;
      m[5334] = 0.000000e+000;
      m[-5334] = 0.000000e+000;
      m[5342] = 1.290893e-012;
      m[-5342] = 1.290893e-012;
      m[5401] = 0.000000e+000;
      m[-5401] = 0.000000e+000;
      m[5403] = 0.000000e+000;
      m[-5403] = 0.000000e+000;
      m[5412] = 1.290893e-012;
      m[-5412] = 1.290893e-012;
      m[5414] = 1.290893e-012;
      m[-5414] = 1.290893e-012;
      m[5422] = 1.290893e-012;
      m[-5422] = 1.290893e-012;
      m[5424] = 1.290893e-012;
      m[-5424] = 1.290893e-012;
      m[5432] = 1.290893e-012;
      m[-5432] = 1.290893e-012;
      m[5434] = 1.290893e-012;
      m[-5434] = 1.290893e-012;
      m[5442] = 1.290893e-012;
      m[-5442] = 1.290893e-012;
      m[5444] = 1.290893e-012;
      m[-5444] = 1.290893e-012;
      m[5503] = 0.000000e+000;
      m[-5503] = 0.000000e+000;
      m[5512] = 1.290893e-012;
      m[-5512] = 1.290893e-012;
      m[5514] = 1.290893e-012;
      m[-5514] = 1.290893e-012;
      m[5522] = 1.290893e-012;
      m[-5522] = 1.290893e-012;
      m[5524] = 1.290893e-012;
      m[-5524] = 1.290893e-012;
      m[5532] = 1.290893e-012;
      m[-5532] = 1.290893e-012;
      m[5534] = 1.290893e-012;
      m[-5534] = 1.290893e-012;
      m[5542] = 1.290893e-012;
      m[-5542] = 1.290893e-012;
      m[5544] = 1.290893e-012;
      m[-5544] = 1.290893e-012;
      m[5554] = 1.290893e-012;
      m[-5554] = 1.290893e-012;
      m[10022] = 0.000000e+000;
      m[10111] = 2.483820e-024;
      m[10113] = 4.635297e-024;
      m[10115] = 2.541360e-024;
      m[10211] = 2.483820e-024;
      m[-10211] = 2.483820e-024;
      m[10213] = 4.635297e-024;
      m[-10213] = 4.635297e-024;
      m[10215] = 2.541360e-024;
      m[-10215] = 2.541360e-024;
      m[9010221] = 1.316424e-023;
      m[10223] = 1.828367e-024;
      m[10225] = 0.000000e+000;
      m[10315] = 3.538775e-024;
      m[-10315] = 3.538775e-024;
      m[10325] = 3.538775e-024;
      m[-10325] = 3.538775e-024;
      m[10331] = 5.265698e-024;
      m[10333] = 0.000000e+000;
      m[10335] = 0.000000e+000;
      m[10443] = 0.000000e+000;
      m[10541] = 0.000000e+000;
      m[-10541] = 0.000000e+000;
      m[10543] = 0.000000e+000;
      m[-10543] = 0.000000e+000;
      m[10551] = 1.000000e+016;
      m[10553] = 0.000000e+000;
      m[10555] = 0.000000e+000;
      m[11112] = 0.000000e+000;
      m[-11112] = 0.000000e+000;
      m[11114] = 2.194041e-024;
      m[-11114] = 2.194041e-024;
      m[11116] = 1.880606e-024;
      m[-11116] = 1.880606e-024;
      m[11212] = 1.880606e-024;
      m[-11212] = 1.880606e-024;
      m[11216] = 0.000000e+000;
      m[-11216] = 0.000000e+000;
      m[12112] = 1.880606e-024;
      m[-12112] = 1.880606e-024;
      m[12114] = 2.194041e-024;
      m[-12114] = 2.194041e-024;
      m[12116] = 5.063171e-024;
      m[-12116] = 5.063171e-024;
      m[12118] = 0.000000e+000;
      m[-12118] = 0.000000e+000;
      m[12122] = 0.000000e+000;
      m[-12122] = 0.000000e+000;
      m[12126] = 1.880606e-024;
      m[-12126] = 1.880606e-024;
      m[12212] = 1.880606e-024;
      m[-12212] = 1.880606e-024;
      m[12214] = 2.194041e-024;
      m[-12214] = 2.194041e-024;
      m[12216] = 5.063171e-024;
      m[-12216] = 5.063171e-024;
      m[12218] = 0.000000e+000;
      m[-12218] = 0.000000e+000;
      m[12222] = 0.000000e+000;
      m[-12222] = 0.000000e+000;
      m[12224] = 2.194041e-024;
      m[-12224] = 2.194041e-024;
      m[12226] = 1.880606e-024;
      m[-12226] = 1.880606e-024;
      m[13112] = 6.582122e-024;
      m[-13112] = 6.582122e-024;
      m[13114] = 1.097020e-023;
      m[-13114] = 1.097020e-023;
      m[13116] = 5.485102e-024;
      m[-13116] = 5.485102e-024;
      m[13122] = 1.316424e-023;
      m[-13122] = 1.316424e-023;
      m[13124] = 1.097020e-023;
      m[-13124] = 1.097020e-023;
      m[13126] = 6.928549e-024;
      m[-13126] = 6.928549e-024;
      m[13212] = 6.582122e-024;
      m[-13212] = 6.582122e-024;
      m[13214] = 1.097020e-023;
      m[-13214] = 1.097020e-023;
      m[13216] = 5.485102e-024;
      m[-13216] = 5.485102e-024;
      m[13222] = 6.582122e-024;
      m[-13222] = 6.582122e-024;
      m[13224] = 1.097020e-023;
      m[-13224] = 1.097020e-023;
      m[13226] = 5.485102e-024;
      m[-13226] = 5.485102e-024;
      m[13314] = 2.742551e-023;
      m[-13314] = 2.742551e-023;
      m[13316] = 0.000000e+000;
      m[-13316] = 0.000000e+000;
      m[13324] = 2.742551e-023;
      m[-13324] = 2.742551e-023;
      m[13326] = 0.000000e+000;
      m[-13326] = 0.000000e+000;
      m[14122] = 1.828367e-022;
      m[-14122] = 1.828367e-022;
      m[14124] = 0.000000e+000;
      m[-14124] = 0.000000e+000;
      m[10221] = 2.194040e-024;
      m[20223] = 2.742551e-023;
      m[20315] = 2.384827e-024;
      m[-20315] = 2.384827e-024;
      m[20325] = 2.384827e-024;
      m[-20325] = 2.384827e-024;
      m[20333] = 1.185968e-023;
      m[20543] = 0.000000e+000;
      m[-20543] = 0.000000e+000;
      m[20553] = 1.000000e+016;
      m[20555] = 0.000000e+000;
      m[21112] = 2.632849e-024;
      m[-21112] = 2.632849e-024;
      m[21114] = 3.291061e-024;
      m[-21114] = 3.291061e-024;
      m[21212] = 2.632849e-024;
      m[-21212] = 2.632849e-024;
      m[21214] = 6.582122e-024;
      m[-21214] = 6.582122e-024;
      m[22112] = 4.388081e-024;
      m[-22112] = 4.388081e-024;
      m[22114] = 3.291061e-024;
      m[-22114] = 3.291061e-024;
      m[22122] = 2.632849e-024;
      m[-22122] = 2.632849e-024;
      m[22124] = 6.582122e-024;
      m[-22124] = 6.582122e-024;
      m[22212] = 4.388081e-024;
      m[-22212] = 4.388081e-024;
      m[22214] = 3.291061e-024;
      m[-22214] = 3.291061e-024;
      m[22222] = 2.632849e-024;
      m[-22222] = 2.632849e-024;
      m[22224] = 3.291061e-024;
      m[-22224] = 3.291061e-024;
      m[23112] = 7.313469e-024;
      m[-23112] = 7.313469e-024;
      m[23114] = 2.991874e-024;
      m[-23114] = 2.991874e-024;
      m[23122] = 4.388081e-024;
      m[-23122] = 4.388081e-024;
      m[23124] = 6.582122e-024;
      m[-23124] = 6.582122e-024;
      m[23126] = 3.291061e-024;
      m[-23126] = 3.291061e-024;
      m[23212] = 7.313469e-024;
      m[-23212] = 7.313469e-024;
      m[23214] = 2.991874e-024;
      m[-23214] = 2.991874e-024;
      m[23222] = 7.313469e-024;
      m[-23222] = 7.313469e-024;
      m[23224] = 2.991874e-024;
      m[-23224] = 2.991874e-024;
      m[23314] = 0.000000e+000;
      m[-23314] = 0.000000e+000;
      m[23324] = 0.000000e+000;
      m[-23324] = 0.000000e+000;
      m[30113] = 2.742551e-024;
      m[30213] = 2.742551e-024;
      m[-30213] = 2.742551e-024;
      m[30223] = 2.991874e-024;
      m[30313] = 2.056913e-024;
      m[-30313] = 2.056913e-024;
      m[30323] = 2.056913e-024;
      m[-30323] = 2.056913e-024;
      m[30343] = 0.000000e+000;
      m[-30343] = 0.000000e+000;
      m[30353] = 0.000000e+000;
      m[-30353] = 0.000000e+000;
      m[30363] = 0.000000e+000;
      m[-30363] = 0.000000e+000;
      m[30411] = 0.000000e+000;
      m[-30411] = 0.000000e+000;
      m[30413] = 0.000000e+000;
      m[-30413] = 0.000000e+000;
      m[30421] = 0.000000e+000;
      m[-30421] = 0.000000e+000;
      m[30423] = 0.000000e+000;
      m[-30423] = 0.000000e+000;
      m[30443] = 2.789035e-023;
      m[30553] = 0.000000e+000;
      m[31114] = 1.880606e-024;
      m[-31114] = 1.880606e-024;
      m[31214] = 4.388081e-024;
      m[-31214] = 4.388081e-024;
      m[32112] = 4.388081e-024;
      m[-32112] = 4.388081e-024;
      m[32114] = 1.880606e-024;
      m[-32114] = 1.880606e-024;
      m[32124] = 4.388081e-024;
      m[-32124] = 4.388081e-024;
      m[32212] = 4.388081e-024;
      m[-32212] = 4.388081e-024;
      m[32214] = 1.880606e-024;
      m[-32214] = 1.880606e-024;
      m[32224] = 1.880606e-024;
      m[-32224] = 1.880606e-024;
      m[33122] = 1.880606e-023;
      m[-33122] = 1.880606e-023;
      m[33314] = 0.000000e+000;
      m[-33314] = 0.000000e+000;
      m[33324] = 0.000000e+000;
      m[-33324] = 0.000000e+000;
      m[41214] = 0.000000e+000;
      m[-41214] = 0.000000e+000;
      m[42112] = 6.582122e-024;
      m[-42112] = 6.582122e-024;
      m[42124] = 0.000000e+000;
      m[-42124] = 0.000000e+000;
      m[42212] = 6.582122e-024;
      m[-42212] = 6.582122e-024;
      m[43122] = 2.194041e-024;
      m[-43122] = 2.194041e-024;
      m[52114] = 0.000000e+000;
      m[-52114] = 0.000000e+000;
      m[52214] = 0.000000e+000;
      m[-52214] = 0.000000e+000;
      m[53122] = 4.388081e-024;
      m[-53122] = 4.388081e-024;
      m[100111] = 1.645531e-024;
      m[100113] = 2.123265e-024;
      m[100211] = 1.645531e-024;
      m[-100211] = 1.645531e-024;
      m[100213] = 2.123265e-024;
      m[-100213] = 2.123265e-024;
      m[100221] = 1.196749e-023;
      m[100223] = 3.871836e-024;
      m[100225] = 0.000000e+000;
      m[100311] = 0.000000e+000;
      m[-100311] = 0.000000e+000;
      m[100313] = 2.837122e-024;
      m[-100313] = 2.837122e-024;
      m[100315] = 0.000000e+000;
      m[-100315] = 0.000000e+000;
      m[100321] = 0.000000e+000;
      m[-100321] = 0.000000e+000;
      m[100323] = 2.837122e-024;
      m[-100323] = 2.837122e-024;
      m[100325] = 0.000000e+000;
      m[-100325] = 0.000000e+000;
      m[100331] = 0.000000e+000;
      m[100333] = 4.388081e-024;
      m[100335] = 3.291061e-024;
      m[100441] = 0.000000e+000;
      m[100551] = 0.000000e+000;
      m[100553] = 1.495937e-020;
      m[100555] = 1.000000e+016;
      m[100557] = 0.000000e+000;
      m[110551] = 1.000000e+016;
      m[110553] = 0.000000e+000;
      m[110555] = 0.000000e+000;
      m[120553] = 1.000000e+016;
      m[120555] = 0.000000e+000;
      m[130553] = 0.000000e+000;
      m[200111] = 3.134344e-024;
      m[200211] = 3.134344e-024;
      m[-200211] = 3.134344e-024;
      m[200551] = 0.000000e+000;
      m[200553] = 2.502708e-020;
      m[200555] = 0.000000e+000;
      m[210551] = 0.000000e+000;
      m[210553] = 0.000000e+000;
      m[220553] = 0.000000e+000;
      m[300553] = 4.701516e-023;
      m[9000221] = 0.000000e+000;
      m[9000443] = 1.265793e-023;
      m[9000553] = 5.983747e-024;
      m[9010443] = 8.438618e-024;
      m[9010553] = 8.331800e-024;
      m[9020221] = 6.038644e-024;
      m[9020443] = 1.530726e-023;
      m[9060225] = 4.388081e-024;
      m[9070225] = 2.056913e-024;
      m[1000001] = 0.000000e+000;
      m[-1000001] = 0.000000e+000;
      m[1000002] = 0.000000e+000;
      m[-1000002] = 0.000000e+000;
      m[1000003] = 0.000000e+000;
      m[-1000003] = 0.000000e+000;
      m[1000004] = 0.000000e+000;
      m[-1000004] = 0.000000e+000;
      m[1000005] = 0.000000e+000;
      m[-1000005] = 0.000000e+000;
      m[1000006] = 0.000000e+000;
      m[-1000006] = 0.000000e+000;
      m[1000011] = 0.000000e+000;
      m[-1000011] = 0.000000e+000;
      m[1000012] = 0.000000e+000;
      m[-1000012] = 0.000000e+000;
      m[1000013] = 0.000000e+000;
      m[-1000013] = 0.000000e+000;
      m[1000014] = 0.000000e+000;
      m[-1000014] = 0.000000e+000;
      m[1000015] = 0.000000e+000;
      m[-1000015] = 0.000000e+000;
      m[1000016] = 0.000000e+000;
      m[-1000016] = 0.000000e+000;
      m[1000021] = 0.000000e+000;
      m[1000022] = 0.000000e+000;
      m[1000023] = 0.000000e+000;
      m[1000024] = 0.000000e+000;
      m[-1000024] = 0.000000e+000;
      m[1000025] = 0.000000e+000;
      m[1000035] = 0.000000e+000;
      m[1000037] = 0.000000e+000;
      m[-1000037] = 0.000000e+000;
      m[1000039] = 0.000000e+000;
      m[2000001] = 0.000000e+000;
      m[-2000001] = 0.000000e+000;
      m[2000002] = 0.000000e+000;
      m[-2000002] = 0.000000e+000;
      m[2000003] = 0.000000e+000;
      m[-2000003] = 0.000000e+000;
      m[2000004] = 0.000000e+000;
      m[-2000004] = 0.000000e+000;
      m[2000005] = 0.000000e+000;
      m[-2000005] = 0.000000e+000;
      m[2000006] = 0.000000e+000;
      m[-2000006] = 0.000000e+000;
      m[2000011] = 0.000000e+000;
      m[-2000011] = 0.000000e+000;
      m[2000012] = 0.000000e+000;
      m[-2000012] = 0.000000e+000;
      m[2000013] = 0.000000e+000;
      m[-2000013] = 0.000000e+000;
      m[2000014] = 0.000000e+000;
      m[-2000014] = 0.000000e+000;
      m[2000015] = 0.000000e+000;
      m[-2000015] = 0.000000e+000;
      m[2000016] = 0.000000e+000;
      m[-2000016] = 0.000000e+000;
      m[3000111] = 0.000000e+000;
      m[3000113] = 0.000000e+000;
      m[3000211] = 0.000000e+000;
      m[-3000211] = 0.000000e+000;
      m[3000213] = 0.000000e+000;
      m[-3000213] = 0.000000e+000;
      m[3000221] = 0.000000e+000;
      m[3000223] = 0.000000e+000;
      m[3000331] = 0.000000e+000;
      m[3100021] = 0.000000e+000;
      m[3100111] = 0.000000e+000;
      m[3100113] = 0.000000e+000;
      m[3200111] = 0.000000e+000;
      m[3200113] = 0.000000e+000;
      m[3300113] = 0.000000e+000;
      m[3400113] = 0.000000e+000;
      m[4000001] = 0.000000e+000;
      m[-4000001] = 0.000000e+000;
      m[4000002] = 0.000000e+000;
      m[-4000002] = 0.000000e+000;
      m[4000011] = 0.000000e+000;
      m[-4000011] = 0.000000e+000;
      m[4000012] = 0.000000e+000;
      m[-4000012] = 0.000000e+000;
      m[5000039] = 0.000000e+000;
      m[9900012] = 0.000000e+000;
      m[9900014] = 0.000000e+000;
      m[9900016] = 0.000000e+000;
      m[9900023] = 0.000000e+000;
      m[9900024] = 0.000000e+000;
      m[-9900024] = 0.000000e+000;
      m[9900041] = 0.000000e+000;
      m[-9900041] = 0.000000e+000;
      m[9900042] = 0.000000e+000;
      m[-9900042] = 0.000000e+000;
      m[1027013000] = 0.000000e+000;
      m[1012006000] = 0.000000e+000;
      m[1063029000] = 0.000000e+000;
      m[1014007000] = 0.000000e+000;
      m[1016008000] = 0.000000e+000;
      m[1028014000] = 0.000000e+000;
      m[1065029000] = 0.000000e+000;
      m[1009004000] = 0.000000e+000;
      m[1019009000] = 0.000000e+000;
      m[1056026000] = 0.000000e+000;
      m[1207082000] = 0.000000e+000;
      m[1208082000] = 0.000000e+000;
      m[1029014000] = 0.000000e+000;
      m[1206082000] = 0.000000e+000;
      m[1054026000] = 0.000000e+000;
      m[1018008000] = 0.000000e+000;
      m[1030014000] = 0.000000e+000;
      m[1057026000] = 0.000000e+000;
      m[1204082000] = 0.000000e+000;
      m[-99000000] = 0.000000e+000;
      m[1028013000] = 0.000000e+000;
      m[1040018000] = 0.000000e+000;
      m[1011005000] = 0.000000e+000;
      m[1012005000] = 0.000000e+000;
      m[1013006000] = 0.000000e+000;
      m[1014006000] = 0.000000e+000;
      m[1052024000] = 0.000000e+000;
      m[1024012000] = 0.000000e+000;
      m[1026012000] = 0.000000e+000;
      m[1027012000] = 0.000000e+000;
      m[1015007000] = 0.000000e+000;
      m[1022010000] = 0.000000e+000;
      m[1058028000] = 0.000000e+000;
      m[1060028000] = 0.000000e+000;
      m[1062028000] = 0.000000e+000;
      m[1064028000] = 0.000000e+000;
      m[1007003000] = 0.000000e+000;
      m[1025012000] = 0.000000e+000;
      m[1053024000] = 0.000000e+000;
      m[1055025000] = 0.000000e+000;
      m[1008004000] = 0.000000e+000;
      m[1010004000] = 0.000000e+000;
      m[1010005000] = 0.000000e+000;
      m[1016007000] = 0.000000e+000;
      m[1017008000] = 0.000000e+000;
      m[1019008000] = 0.000000e+000;
      m[1023010000] = 0.000000e+000;
      m[1024011000] = 0.000000e+000;
      m[1031015000] = 0.000000e+000;
      m[1039017000] = 0.000000e+000;
      m[1040017000] = 0.000000e+000;
      m[1036018000] = 0.000000e+000;
      m[1050024000] = 0.000000e+000;
      m[1054024000] = 0.000000e+000;
      m[1059026000] = 0.000000e+000;
      m[1061028000] = 0.000000e+000;
      m[1063028000] = 0.000000e+000;
      m[1092042000] = 0.000000e+000;
      m[1095042000] = 0.000000e+000;
      m[1096042000] = 0.000000e+000;
      m[1097042000] = 0.000000e+000;
      m[1098042000] = 0.000000e+000;
      m[1100042000] = 0.000000e+000;
      m[1108046000] = 0.000000e+000;

      // Added by hand:
      m[9902210] = 0.000000e+000; //diffractive p-state -> assume no lifetime
      return true;
    }


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_mult_total;  // full kinematic range
    Histo1DPtr _h_mult_eta[5]; // in eta bins
    Histo1DPtr _h_mult_pt[5];  // in pT bins
    Histo1DPtr _h_dndeta;       // density dn/deta
    Histo1DPtr _h_dndpt;    // density dn/dpT
    //@}


    /// @name Private variables
    double _p_min;
    double _pt_min;
    double _eta_min;
    double _eta_max;
    double _maxlft;

    /// Count selected events
    CounterPtr _sumW;

    map<int, double> _partLftMap; // Map <PDGID, PDGLIFETIME>

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCB_2014_I1281685);

}
