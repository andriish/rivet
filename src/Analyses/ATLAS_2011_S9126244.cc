// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  struct ATLAS_2011_S9126244_Plots {

    int selectionType; // The third value in the axis code d**-**-**
    std::string intermediateHistName;

    // Gap fraction vs DeltaY plot setup
    int m_gapFractionDeltaYHistIndex;
    std::vector<double> m_gapFractionDeltaYSlices;
    BinnedHistogram<double> _h_gapVsDeltaYVeto;
    BinnedHistogram<double> _h_gapVsDeltaYInc;

    // Gap fraction vs ptBar plot setup
    int m_gapFractionPtBarHistIndex;
    std::vector<double> m_gapFractionPtBarSlices;
    BinnedHistogram<double> _h_gapVsPtBarVeto;
    BinnedHistogram<double> _h_gapVsPtBarInc;

    // Gap fraction vs Q0 plot setup
    int m_gapFractionQ0HistIndex;
    std::vector<double> m_gapFractionQ0SlicesPtBar;
    std::vector<double> m_gapFractionQ0SlicesDeltaY;
    std::vector<Histo1DPtr> _h_vetoPt;
    std::vector<Scatter2DPtr> _d_vetoPtGapFraction;
    std::vector<double> _h_vetoPtTotalSum;

    // Average njet vs DeltaY setup
    int m_avgNJetDeltaYHistIndex;
    std::vector<double> m_avgNJetDeltaYSlices;
    std::vector<Profile1DPtr> _p_avgJetVsDeltaY;

    // Average njet vs PptBar setup
    int m_avgNJetPtBarHistIndex;
    std::vector<double> m_avgNJetPtBarSlices;
    std::vector<Profile1DPtr> _p_avgJetVsPtBar;
  };



  /// ATLAS dijet production with central jet veto
  /// @todo Need work to make sure that temp histos are removed
  class ATLAS_2011_S9126244 : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_S9126244()
      : Analysis("ATLAS_2011_S9126244")
    {    }


  public:

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialize the lone projection required
      addProjection(FastJets(FinalState(), FastJets::ANTIKT, 0.6), "AntiKtJets06");

      // Make Q0 bins 0->20 then 20 to 195.0 in steps of 5
      m_q0BinEdges += 0.0;
      for (size_t x = 0; x < 36; x++) {
        m_q0BinEdges += 20.0 + x*5.0;
      }

      // Initialize plots for each selection type
      m_selectionPlots[0].intermediateHistName = "highestPt";
      m_selectionPlots[0].selectionType = 1;
      m_selectionPlots[0].m_gapFractionDeltaYHistIndex = 6;
      m_selectionPlots[0].m_gapFractionPtBarHistIndex = 1;
      m_selectionPlots[0].m_gapFractionQ0HistIndex = 13;
      m_selectionPlots[0].m_avgNJetDeltaYHistIndex = 37;
      m_selectionPlots[0].m_avgNJetPtBarHistIndex = 26;
      m_selectionPlots[0].m_gapFractionDeltaYSlices += 70.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0;
      m_selectionPlots[0].m_gapFractionPtBarSlices += 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
      m_selectionPlots[0].m_gapFractionQ0SlicesPtBar += 70.0, 90.0, 120.0, 150.0, 210.0, 240.0;
      m_selectionPlots[0].m_gapFractionQ0SlicesDeltaY += 2.0, 3.0, 4.0, 5.0;
      m_selectionPlots[0].m_avgNJetPtBarSlices += 1.0, 2.0, 3.0, 4.0, 5.0;
      m_selectionPlots[0].m_avgNJetDeltaYSlices += 70.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0;
      initializePlots(m_selectionPlots[0]);

      m_selectionPlots[1].intermediateHistName = "forwardBackward";
      m_selectionPlots[1].selectionType = 2;
      m_selectionPlots[1].m_gapFractionDeltaYHistIndex = 6;
      m_selectionPlots[1].m_gapFractionPtBarHistIndex = 1;
      m_selectionPlots[1].m_gapFractionQ0HistIndex = 13;
      m_selectionPlots[1].m_avgNJetDeltaYHistIndex = 37;
      m_selectionPlots[1].m_avgNJetPtBarHistIndex = 26;
      m_selectionPlots[1].m_gapFractionDeltaYSlices += 70.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0;
      m_selectionPlots[1].m_gapFractionPtBarSlices += 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
      m_selectionPlots[1].m_gapFractionQ0SlicesPtBar += 70.0, 90.0, 120.0, 150.0, 210.0, 240.0;
      m_selectionPlots[1].m_gapFractionQ0SlicesDeltaY += 2.0, 3.0, 4.0, 5.0;
      m_selectionPlots[1].m_avgNJetPtBarSlices += 1.0, 2.0, 3.0, 4.0, 5.0;
      m_selectionPlots[1].m_avgNJetDeltaYSlices += 70.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0;
      initializePlots(m_selectionPlots[1]);

      m_selectionPlots[2].intermediateHistName = "forwardBackward_PtBarVeto";
      m_selectionPlots[2].selectionType = 1;
      m_selectionPlots[2].m_gapFractionDeltaYHistIndex = 19;
      m_selectionPlots[2].m_avgNJetDeltaYHistIndex = 30;
      m_selectionPlots[2].m_gapFractionDeltaYSlices += 70.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0;
      m_selectionPlots[2].m_avgNJetDeltaYSlices += 70.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0;
      initializePlots(m_selectionPlots[2]);
    }


    /// @todo Replace with unregistered temp histos, maybe binned by YODA -- or use automatic discarding of TMP histos.
    void initializePlots(ATLAS_2011_S9126244_Plots& plots) {

      // Gap fraction vs DeltaY
      if (!plots.m_gapFractionDeltaYSlices.empty()) {
      for (size_t x = 0; x < plots.m_gapFractionDeltaYSlices.size()-1; x++) {
        const string vetoHistName = "gapDeltaYVeto_" + plots.intermediateHistName + "_" + lexical_cast<string>(x);
        const string inclusiveHistName = "gapDeltaYInclusive_" + plots.intermediateHistName + "_" + lexical_cast<string>(x);
        plots._h_gapVsDeltaYVeto.addHistogram(plots.m_gapFractionDeltaYSlices[x], plots.m_gapFractionDeltaYSlices[x+1],
                                              bookHisto1D(plots.m_gapFractionDeltaYHistIndex+x, 1, plots.selectionType, vetoHistName));
        plots._h_gapVsDeltaYInc.addHistogram(plots.m_gapFractionDeltaYSlices[x], plots.m_gapFractionDeltaYSlices[x+1],
                                             bookHisto1D(plots.m_gapFractionDeltaYHistIndex+x, 1, plots.selectionType, inclusiveHistName));
      }
      }

      // Average njet vs DeltaY
      if (!plots.m_avgNJetDeltaYSlices.empty()) {
      for (size_t x = 0; x < plots.m_avgNJetDeltaYSlices.size()-1; x++) {
        plots._p_avgJetVsDeltaY += bookProfile1D(plots.m_avgNJetDeltaYHistIndex+x, 1, plots.selectionType);
      }
      }

      // Gap fraction vs PtBar
      if (!plots.m_gapFractionPtBarSlices.empty()) {
      for (size_t x = 0; x < plots.m_gapFractionPtBarSlices.size()-1; x++) {
        const string vetoHistName = "gapPtBarVeto_" + plots.intermediateHistName + "_" + lexical_cast<string>(x);
        const string inclusiveHistName = "gapPtBarInclusive_" + plots.intermediateHistName + "_" + lexical_cast<string>(x);
        plots._h_gapVsPtBarVeto.addHistogram(plots.m_gapFractionPtBarSlices[x], plots.m_gapFractionPtBarSlices[x+1],
                                             bookHisto1D(plots.m_gapFractionPtBarHistIndex+x, 1, plots.selectionType, vetoHistName));
        plots._h_gapVsPtBarInc.addHistogram(plots.m_gapFractionPtBarSlices[x], plots.m_gapFractionPtBarSlices[x+1],
                                            bookHisto1D(plots.m_gapFractionPtBarHistIndex+x, 1, plots.selectionType, inclusiveHistName));
      }
      }

      // Average njet vs PtBar
      if (!plots.m_avgNJetPtBarSlices.empty()) {
      for (size_t x=0; x<plots.m_avgNJetPtBarSlices.size()-1; x++) {
        plots._p_avgJetVsPtBar += bookProfile1D(plots.m_avgNJetPtBarHistIndex+x, 1, plots.selectionType);
      }
      }

      // Gap fraction vs Q0
      int q0PlotCount = 0;
      for (size_t x = 0; x < plots.m_gapFractionQ0SlicesPtBar.size()/2; x++) {
        for (size_t y = 0; y < plots.m_gapFractionQ0SlicesDeltaY.size()/2; y++) {
          const string vetoPtHistName = "vetoPt_" + plots.intermediateHistName + "_" + lexical_cast<string>(q0PlotCount);
          const string vetoPtGapDataPointName = "gapQ0GapFractionDataPoints_" + plots.intermediateHistName + "_" + lexical_cast<string>(q0PlotCount);
          plots._h_vetoPt += bookHisto1D(vetoPtHistName, m_q0BinEdges);
          plots._d_vetoPtGapFraction += bookScatter2D(plots.m_gapFractionQ0HistIndex+q0PlotCount, 1, plots.selectionType);
          plots._h_vetoPtTotalSum += 0.0;
          q0PlotCount++;
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get minimal list of jets needed to be considered
      double minimumJetPtBar = 50.0*GeV; // of interval defining jets

      vector<FourMomentum> acceptJets;
      foreach (const Jet& jet, applyProjection<FastJets>(event, "AntiKtJets06").jetsByPt(20.0*GeV)) {
        if (fabs(jet.momentum().rapidity()) < 4.4) {
          acceptJets.push_back(jet.momentum());
        }
      }

      // If we can't form an interval, drop out of the analysis early
      if (acceptJets.size() < 2) vetoEvent;

      // Analyze leading jet case
      if (acceptJets[0].pT() + acceptJets[1].pT() > 2*minimumJetPtBar) {
        analyzeJets(acceptJets, m_selectionPlots[0], weight, 20.0*GeV);
      }

      // Find the most forward-backward jets
      size_t minRapidityJet = 0, maxRapidityJet = 0;
      for (size_t j = 1; j < acceptJets.size(); j++) {
        if (acceptJets[j].rapidity() > acceptJets[maxRapidityJet].rapidity()) maxRapidityJet = j;
        if (acceptJets[j].rapidity() < acceptJets[minRapidityJet].rapidity()) minRapidityJet = j;
      }

      // Make a container of jet momenta with the extreme f/b jets at the front
      vector<FourMomentum> fwdBkwdJets;
      fwdBkwdJets.push_back(acceptJets[maxRapidityJet]);
      fwdBkwdJets.push_back(acceptJets[minRapidityJet]);
      for (size_t j = 0; j < acceptJets.size(); j++) {
        if (j == minRapidityJet || j == maxRapidityJet) continue;
        fwdBkwdJets.push_back(acceptJets[j]);
      }

      if (fwdBkwdJets[0].pT() + fwdBkwdJets[1].pT() > 2*minimumJetPtBar) {
        // Use most forward/backward jets in rapidity to define the interval
        analyzeJets(fwdBkwdJets, m_selectionPlots[1], weight, 20.0*GeV);
        // As before but now using PtBar of interval to define veto threshold
        analyzeJets(fwdBkwdJets, m_selectionPlots[2], weight, (fwdBkwdJets[0].pT()+fwdBkwdJets[1].pT())/2.0);
      }
    }


    /// Fill plots!
    void analyzeJets(vector<FourMomentum>& jets, ATLAS_2011_S9126244_Plots& plots,
                     const double weight, double vetoPtThreshold) {

      // Calculate the interval size, ptBar and veto Pt (if any)
      const double intervalSize = fabs(jets[0].rapidity()-jets[1].rapidity());
      const double ptBar = (jets[0].pT()+jets[1].pT())/2.0;

      const double minY = min(jets[0].rapidity(), jets[1].rapidity());
      const double maxY = max(jets[0].rapidity(), jets[1].rapidity());

      double vetoPt = 0.0*GeV;
      for (size_t j = 2; j < jets.size(); j++) {
        if (inRange(jets[j].rapidity(), minY, maxY)) vetoPt = max(jets[j].pT(), vetoPt);
      }

      // Fill the gap fraction vs delta Y histograms
      plots._h_gapVsDeltaYInc.fill(ptBar/GeV, intervalSize, weight);
      if (vetoPt < vetoPtThreshold) {
        plots._h_gapVsDeltaYVeto.fill(ptBar/GeV, intervalSize, weight);
      }

      // Fill the gap fraction vs ptBar histograms
      plots._h_gapVsPtBarInc.fill(intervalSize, ptBar/GeV,  weight);
      if (vetoPt < vetoPtThreshold) {
        plots._h_gapVsPtBarVeto.fill(intervalSize, ptBar/GeV, weight);
      }

      // Count the number of veto jets present
      int vetoJetsCount = 0;
      for (size_t j = 2; j < jets.size(); j++) {
        if (inRange(jets[j].rapidity(), minY, maxY) && jets[j].pT() > vetoPtThreshold) {
          vetoJetsCount += 1;
        }
      }

      // Fill the avg NJet, deltaY slices
      if (!plots.m_avgNJetPtBarSlices.empty()) {
      for (size_t i = 0; i < plots.m_avgNJetPtBarSlices.size()-1; i++) {
        if (inRange(intervalSize, plots.m_avgNJetPtBarSlices[i], plots.m_avgNJetPtBarSlices[i+1])) {
          plots._p_avgJetVsPtBar[i]->fill(ptBar/GeV, vetoJetsCount, weight);
        }
      }
      }

      // Fill the avg NJet, ptBar slices
      if (!plots.m_avgNJetDeltaYSlices.empty()) {
      for (size_t i = 0; i < plots.m_avgNJetDeltaYSlices.size()-1; i++) {
        if (inRange(ptBar/GeV, plots.m_avgNJetDeltaYSlices[i], plots.m_avgNJetDeltaYSlices[i+1])) {
          plots._p_avgJetVsDeltaY[i]->fill(intervalSize, vetoJetsCount, weight);
        }
      }
      }

      // Fill the veto pt plots
      int q0PlotCount = 0;
      for (size_t x = 0; x < plots.m_gapFractionQ0SlicesPtBar.size()/2; x++) {
        for (size_t y = 0; y < plots.m_gapFractionQ0SlicesDeltaY.size()/2; y++) {
          // Check if it should be filled
          if ( ptBar/GeV < plots.m_gapFractionQ0SlicesPtBar[x*2] ||
               ptBar/GeV >= plots.m_gapFractionQ0SlicesPtBar[x*2+1] ) {
            q0PlotCount++;
            continue;
          }

          if ( intervalSize < plots.m_gapFractionQ0SlicesDeltaY[y*2] ||
               intervalSize >= plots.m_gapFractionQ0SlicesDeltaY[y*2+1] ) {
            q0PlotCount++;
            continue;
          }

          plots._h_vetoPt[q0PlotCount]->fill(vetoPt, weight);
          plots._h_vetoPtTotalSum[q0PlotCount] += weight;

          q0PlotCount++;
        }
      }
    }


    /// Derive final distributions for each selection
    void finalize() {
      foreach (const ATLAS_2011_S9126244_Plots& plots, m_selectionPlots) {

        /// @todo Clean up temp histos -- requires restructuring the temp histo struct

        for (size_t x = 0; x < plots._h_gapVsDeltaYVeto.getHistograms().size(); x++) {
          divide(plots._h_gapVsDeltaYVeto.getHistograms()[x], plots._h_gapVsDeltaYInc.getHistograms()[x],
                 bookScatter2D(plots.m_gapFractionDeltaYHistIndex+x, 1, plots.selectionType));
        }
        for (size_t x = 0; x < plots._h_gapVsPtBarVeto.getHistograms().size(); x++) {
          divide(plots._h_gapVsPtBarVeto.getHistograms()[x], plots._h_gapVsPtBarInc.getHistograms()[x],
                 bookScatter2D(plots.m_gapFractionPtBarHistIndex+x, 1, plots.selectionType));
        }

        for (size_t h = 0; h < plots._d_vetoPtGapFraction.size(); h++) {
          size_t numbins = refData(plots.m_gapFractionQ0HistIndex+h, 1, plots.selectionType).numPoints();
          finalizeQ0GapFraction(plots._h_vetoPtTotalSum[h],
                                plots._d_vetoPtGapFraction[h],
                                plots._h_vetoPt[h], numbins);
         }
      }
    }


    /// Convert the differential histograms to an integral histo and assign binomial errors as a efficiency
    /// @todo Should be convertible to a YODA ~one-liner
    void finalizeQ0GapFraction(double totalWeightSum, Scatter2DPtr gapFractionDP, Histo1DPtr vetoPtHist, size_t numbins) {
      double vetoPtWeightSum = 0.0;
      for (size_t i = 0; i < numbins; i++) {
        vetoPtWeightSum += vetoPtHist->bin(i).sumW();

        // Calculate the efficiency uncertainty
        const double eff = (totalWeightSum != 0) ? vetoPtWeightSum/totalWeightSum : 0;
        const double effErr = (totalWeightSum != 0) ? sqrt( eff*(1.0-eff)/totalWeightSum ) : 0;
        gapFractionDP->point(i).setY(eff, effErr);
      }
      /// @todo Remove vetoPtHist somehow
    }


  private:

    // Only need to define the q0 binning here
    vector<double> m_q0BinEdges;

    // Struct containing complete set of plots, times 3
    ATLAS_2011_S9126244_Plots m_selectionPlots[3];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_S9126244);

}
