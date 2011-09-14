// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  struct ATLAS_2011_S9126244_Plots {

    int selectionType; // The third value in the axis code d**-**-**
    std::string intermediateHistName;

    // Gap Fraction Vs Delta Y Plot Setup
    int m_gapFractionDeltaYHistIndex;
    std::vector<double> m_gapFractionDeltaYSlices;
    BinnedHistogram<double> _h_gapVsDeltaYVeto;
    BinnedHistogram<double> _h_gapVsDeltaYInc;

    // Gap Fraction Vs PtBar Plot Setup
    int m_gapFractionPtBarHistIndex;
    std::vector<double> m_gapFractionPtBarSlices;
    BinnedHistogram<double> _h_gapVsPtBarVeto;
    BinnedHistogram<double> _h_gapVsPtBarInc;

    // Gap Fraction Vs Q0 Plot Setup
    int m_gapFractionQ0HistIndex;
    std::vector<double> m_gapFractionQ0SlicesPtBar;
    std::vector<double> m_gapFractionQ0SlicesDeltaY;
    std::vector<AIDA::IHistogram1D*> _h_vetoPt;
    std::vector<AIDA::IDataPointSet*> _d_vetoPtGapFraction;
    std::vector<double> _h_vetoPtTotalSum;

    // Average NJet Vs DeltaY Setup
    int m_avgNJetDeltaYHistIndex;
    std::vector<double> m_avgNJetDeltaYSlices;
    std::vector<AIDA::IProfile1D*> _p_avgJetVsDeltaY;

    // Average NJet Vs PtBar Setup
    int m_avgNJetPtBarHistIndex;
    std::vector<double> m_avgNJetPtBarSlices;
    std::vector<AIDA::IProfile1D*> _p_avgJetVsPtBar;
  };


  class ATLAS_2011_S9126244 : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_S9126244()
      : Analysis("ATLAS_2011_S9126244")
    {
      setNeedsCrossSection(false);
    }


  public:

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialize the lone projection required
      addProjection(FastJets(FinalState(), FastJets::ANTIKT, 0.6), "AntiKtJets06");

      // Make Q0 bins 0->20 then 20 to 195.0 in steps of 5
      m_q0BinEdges += 0.0;
      for (unsigned int x=0; x<36; x++){
        m_q0BinEdges += 20.0 + x*5.0 ;
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


    void initializePlots(ATLAS_2011_S9126244_Plots& plots) {

      // Gap Fraction Vs DeltaY
      for (int x = 0; x < ((int)plots.m_gapFractionDeltaYSlices.size()-1); x++) {
        std::stringstream vetoHistName;
        std::stringstream inclusiveHistName;
        const BinEdges deltaYEdges = binEdges(plots.m_gapFractionDeltaYHistIndex+x, 1, plots.selectionType);

        vetoHistName << "gapDeltaYVeto_" << plots.intermediateHistName << "_" << x;
        inclusiveHistName << "gapDeltaYInclusive_" << plots.intermediateHistName << "_" << x;

        plots._h_gapVsDeltaYVeto.addHistogram(plots.m_gapFractionDeltaYSlices[x], plots.m_gapFractionDeltaYSlices[x+1], bookHistogram1D(vetoHistName.str(), deltaYEdges));
        plots._h_gapVsDeltaYInc.addHistogram(plots.m_gapFractionDeltaYSlices[x], plots.m_gapFractionDeltaYSlices[x+1], bookHistogram1D(inclusiveHistName.str(), deltaYEdges));
      }

      // Average NJet Vs DeltaY
      for (int x = 0; x < ((int)plots.m_avgNJetDeltaYSlices.size()-1); x++) {
        plots._p_avgJetVsDeltaY += bookProfile1D(plots.m_avgNJetDeltaYHistIndex+x, 1, plots.selectionType);
      }

      // Gap Fraction Vs PtBar
      for (int x = 0; x < ((int)plots.m_gapFractionPtBarSlices.size()-1); x++) {

        std::stringstream vetoHistName;
        std::stringstream inclusiveHistName;
        const BinEdges ptBarEdges = binEdges(plots.m_gapFractionPtBarHistIndex+x, 1, plots.selectionType);

        vetoHistName << "gapPtBarVeto_" << plots.intermediateHistName << "_" << x;
        inclusiveHistName << "gapPtBarInclusive_" << plots.intermediateHistName << "_" << x;

        plots._h_gapVsPtBarVeto.addHistogram(plots.m_gapFractionPtBarSlices[x], plots.m_gapFractionPtBarSlices[x+1], bookHistogram1D(vetoHistName.str(), ptBarEdges));
        plots._h_gapVsPtBarInc.addHistogram(plots.m_gapFractionPtBarSlices[x], plots.m_gapFractionPtBarSlices[x+1], bookHistogram1D(inclusiveHistName.str(), ptBarEdges));
      }

      // Average NJet Vs PtBar
      for (int x = 0; x < ((int)plots.m_avgNJetPtBarSlices.size()-1); x++) {
        plots._p_avgJetVsPtBar += bookProfile1D(plots.m_avgNJetPtBarHistIndex+x, 1, plots.selectionType);
      }

      // Gap fraction Vs Q0
      int q0PlotCount = 0;
      for (int x = 0; x < ((int)plots.m_gapFractionQ0SlicesPtBar.size()/2); x++) {
        for (int y = 0; y < ((int)plots.m_gapFractionQ0SlicesDeltaY.size()/2); y++) {
          std::stringstream vetoPtHistName;
          std::stringstream vetoPtGapDataPointName;

          vetoPtHistName << "vetoPt_" << plots.intermediateHistName << "_" << q0PlotCount;
          vetoPtGapDataPointName << "gapQ0GapFractionDataPoints_" << plots.intermediateHistName << "_" << q0PlotCount;

          plots._h_vetoPt += bookHistogram1D(vetoPtHistName.str(),
                                             m_q0BinEdges);
          plots._d_vetoPtGapFraction += bookDataPointSet(plots.m_gapFractionQ0HistIndex+q0PlotCount, 1, plots.selectionType);
          plots._h_vetoPtTotalSum += 0.0;
          q0PlotCount++;
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get minimal list of jets needed to be considered
      double minimumJetPtBar = 50.0*GeV; //of interval defining jets

      vector<FourMomentum> acceptJets;
      foreach (const Jet& jet, applyProjection<FastJets>(event, "AntiKtJets06").jetsByPt(20.0*GeV)) {
        if (fabs(jet.momentum().rapidity()) < 4.4) {
          acceptJets.push_back(jet.momentum());
        }
      }

      // If can't form an interval drop out of the analysis early
      if (acceptJets.size() < 2) {
        return;
      }

      // Analyze leading jet case
      if ((acceptJets[0].pT() + acceptJets[1].pT())/2.0 > minimumJetPtBar) {
        analyzeJets(acceptJets, m_selectionPlots[0], weight, 20.0*GeV);
      }

      // Re-order jets to have forward backward selection
      unsigned int minRapidityJet = 0;
      unsigned int maxRapidityJet = 0;

      for (size_t j=1; j<acceptJets.size(); j++) {
        if (acceptJets[j].rapidity() >
            acceptJets[maxRapidityJet].rapidity()) {
          maxRapidityJet=j;
        }

        if (acceptJets[j].rapidity() <
            acceptJets[minRapidityJet].rapidity()) {
          minRapidityJet=j;
        }
      }

      vector<FourMomentum> fwdBkwdJets;
      fwdBkwdJets.push_back(acceptJets[maxRapidityJet]);
      fwdBkwdJets.push_back(acceptJets[minRapidityJet]);

      for (size_t j=0; j<acceptJets.size(); j++) {
        if (j==minRapidityJet or j==maxRapidityJet){
          continue;
        }
        fwdBkwdJets.push_back(acceptJets[j]);
      }

      if ((fwdBkwdJets[0].pT() + fwdBkwdJets[1].pT())/2.0 > minimumJetPtBar) {
        //Use most forward/backward jets in rapidity to define the interval
        analyzeJets(fwdBkwdJets, m_selectionPlots[1], weight,
                    20.0*GeV);

        //As before but now using PtBar of interval to define veto threshold
        analyzeJets(fwdBkwdJets, m_selectionPlots[2], weight,
                    (fwdBkwdJets[0].pT()+fwdBkwdJets[1].pT())/2.0);
      }
    }


    // Fill plots!
    void analyzeJets(vector<FourMomentum>& jets, ATLAS_2011_S9126244_Plots& plots,
                     const double weight, double vetoPtThreshold){
      // Calculate the interval size, ptBar and veto Pt (if any)
      double intervalSize;
      double ptBar;
      double vetoPt = 0.0*GeV;

      intervalSize = fabs(jets[0].rapidity()-jets[1].rapidity());
      ptBar = (jets[0].pT()+jets[1].pT())/2.0;

      double minY;
      double maxY;
      if (jets[0].rapidity() > jets[1].rapidity()) {
        minY = jets[1].rapidity();
        maxY = jets[0].rapidity();
      } else {
        minY = jets[0].rapidity();
        maxY = jets[1].rapidity();
      }

      for (size_t j=2; j<jets.size(); j++){
        if (jets[j].rapidity() > minY &&
            jets[j].rapidity() < maxY &&
            jets[j].pT() > vetoPt){
          vetoPt = jets[j].pT();
        }
      }

      // Fill the gap fraction vs delta Y histograms
      plots._h_gapVsDeltaYInc.fill(ptBar/GeV, intervalSize, weight);
      if (vetoPt < vetoPtThreshold) {
        plots._h_gapVsDeltaYVeto.fill(ptBar/GeV, intervalSize, weight);
      }

      // Fill the gap fraction vs pt Bar histograms
      plots._h_gapVsPtBarInc.fill(intervalSize, ptBar/GeV,  weight);
      if (vetoPt < vetoPtThreshold) {
        plots._h_gapVsPtBarVeto.fill(intervalSize, ptBar/GeV, weight);
      }

      // Count the number of veto jets present
      int vetoJetsCount=0;
      for (size_t j=2; j<jets.size(); j++){
        if (jets[j].rapidity() > minY &&
            jets[j].rapidity() < maxY &&
            jets[j].pT() > vetoPtThreshold){
          vetoJetsCount += 1;
        }
      }

      // Fill the avg NJet, deltaY slices
      for (int i=0; i<(int)plots.m_avgNJetPtBarSlices.size()-1; i++) {
        if ( intervalSize >= plots.m_avgNJetPtBarSlices[i] &&
             intervalSize < plots.m_avgNJetPtBarSlices[i+1]) {
         plots._p_avgJetVsPtBar[i]->fill(ptBar/GeV, vetoJetsCount, weight);
        }
      }

      // Fill the avg NJet, ptBar slices
      for (int i=0; i<(int)plots.m_avgNJetDeltaYSlices.size()-1; i++) {
        if ( ptBar/GeV >= plots.m_avgNJetDeltaYSlices[i] &&
             ptBar/GeV < plots.m_avgNJetDeltaYSlices[i+1] ) {
          plots._p_avgJetVsDeltaY[i]->fill(intervalSize, vetoJetsCount, weight);
        }
      }

      // Fill the veto pt plots
      int q0PlotCount = 0;
      for (int x=0; x<((int)plots.m_gapFractionQ0SlicesPtBar.size()/2); x++) {
        for (int y=0; y<((int)plots.m_gapFractionQ0SlicesDeltaY.size()/2); y++) {
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
        // Calculate the gap fraction for each slice
        for (size_t x = 0; x < plots._h_gapVsDeltaYVeto.getHistograms().size(); x++) {
          histogramFactory().divide(histoPath(makeAxisCode(plots.m_gapFractionDeltaYHistIndex+x, 1, plots.selectionType)),
                                    *(plots._h_gapVsDeltaYVeto.getHistograms()[x]),
                                    *(plots._h_gapVsDeltaYInc.getHistograms()[x]));
          histogramFactory().destroy(plots._h_gapVsDeltaYVeto.getHistograms()[x]);
          histogramFactory().destroy(plots._h_gapVsDeltaYInc.getHistograms()[x]);
        }

        for (size_t x = 0; x < plots._h_gapVsPtBarVeto.getHistograms().size(); x++) {
          histogramFactory().divide(histoPath(makeAxisCode(plots.m_gapFractionPtBarHistIndex+x, 1, plots.selectionType)),
                                    *(plots._h_gapVsPtBarVeto.getHistograms()[x]),
                                    *(plots._h_gapVsPtBarInc.getHistograms()[x]));
          histogramFactory().destroy(plots._h_gapVsPtBarVeto.getHistograms()[x]);
          histogramFactory().destroy(plots._h_gapVsPtBarInc.getHistograms()[x]);
        }

        for (size_t h = 0; h < plots._d_vetoPtGapFraction.size(); h++) {
          // Get the number of bins needed for this slice
          const BinEdges q0Edges = binEdges(plots.m_gapFractionQ0HistIndex+h, 1, plots.selectionType);
          finalizeQ0GapFraction(plots._h_vetoPtTotalSum[h],
                                plots._d_vetoPtGapFraction[h],
                                plots._h_vetoPt[h],
                                q0Edges.size());
        }
      }
    }


    void finalizeQ0GapFraction(double totalWeightSum,
                               AIDA::IDataPointSet* gapFractionDP,
                               AIDA::IHistogram1D* vetoPtHist,
                               int binNumber) {
      double vetoPtWeightSum = 0.0;
      for (int x = 0; x < binNumber-1; x++) {
        vetoPtWeightSum += vetoPtHist->binHeight(x);

        // Alternatively try saving as data points
        IDataPoint* currentPoint = gapFractionDP->point(x);
        IMeasurement* xCoord = currentPoint->coordinate(0);
        IMeasurement* yCoord = currentPoint->coordinate(1);

        // Calculate the efficiency uncertainty
        double efficiency = vetoPtWeightSum/totalWeightSum;
        double efficiencyError = std::sqrt(efficiency*(1.0-efficiency)/totalWeightSum);
        if (totalWeightSum==0.) efficiency = efficiencyError = 0.;

        xCoord->setValue(m_q0BinEdges[x+1]);
        xCoord->setErrorPlus(2.5);
        xCoord->setErrorMinus(2.5);
        yCoord->setValue(efficiency);
        yCoord->setErrorPlus(efficiencyError);
        yCoord->setErrorMinus(efficiencyError);
      }
      histogramFactory().destroy(vetoPtHist);
    }


  private:

    // Only need to define the q0 binning here
    std::vector<double> m_q0BinEdges;


  private:

    // Structure containing complete set of plots
    ATLAS_2011_S9126244_Plots m_selectionPlots[3];

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_S9126244);

}
