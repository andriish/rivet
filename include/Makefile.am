SUBDIRS = Rivet

#dist_noinst_HEADERS = AIDA/AIDA.h AIDA/IAnalysisFactory.h \
#                      AIDA/ITree.h AIDA/IHistogramFactory.h 

nobase_include_HEADERS =

if NO_INSTALLED_AIDA
nobase_include_HEADERS += \
  AIDA/IFitParameterSettings.h    AIDA/IHistogram.h               AIDA/IDataPointSetFactory.h     AIDA/ICloud3D.h  \
  AIDA/IFunctionCatalog.h         AIDA/IDataPointSet.h            AIDA/IPlotterFactory.h          AIDA/IPlotterRegion.h  \
  AIDA/ITupleFactory.h            AIDA/ITree.h                    AIDA/IDataPoint.h               AIDA/IEvaluator.h  \
  AIDA/IHistogram2D.h             AIDA/IProfile.h                 AIDA/IManagedObject.h           AIDA/IMarkerStyle.h  \
  AIDA/IFunction.h                AIDA/ITitleStyle.h              AIDA/IHistogram3D.h             AIDA/IAxis.h  \
  AIDA/AIDA.h                     AIDA/IFunctionFactory.h         AIDA/IFilter.h                  AIDA/IFillStyle.h  \
  AIDA/IAxisStyle.h               AIDA/IBrushStyle.h              AIDA/IHistogramFactory.h        AIDA/IRangeSet.h  \
  AIDA/IProfile2D.h               AIDA/ITextStyle.h               AIDA/IMeasurement.h             AIDA/IFitFactory.h  \
  AIDA/IPlotterLayout.h           AIDA/ICloud.h                   AIDA/IProfile1D.h               AIDA/IAnnotation.h  \
  AIDA/IBaseHistogram.h           AIDA/IModelFunction.h           AIDA/IAnalysisFactory.h         AIDA/IPlotterStyle.h  \
  AIDA/IPlotter.h                 AIDA/ILineStyle.h               AIDA/ITreeFactory.h             AIDA/IFitData.h  \
  AIDA/IDataStyle.h               AIDA/ITupleEntry.h              AIDA/ITuple.h                   AIDA/IInfoStyle.h  \
  AIDA/IBaseStyle.h               AIDA/IConstants.h               AIDA/ICloud1D.h                 AIDA/IInfo.h  \
  AIDA/IHistogram1D.h             AIDA/ICloud2D.h                 AIDA/IFitResult.h               AIDA/IFitter.h
endif

if NO_INSTALLED_LWH
nobase_include_HEADERS += \
  LWH/AIAnalysisFactory.h         LWH/AIManagedObject.h           LWH/Axis.h                      LWH/VariAxis.h  \
  LWH/AIAxis.h                    LWH/AIDataPoint.h               LWH/Tree.h                      LWH/DataPoint.h  \
  LWH/AIHistogram1D.h             LWH/AIDataPointSet.h            LWH/TreeFactory.h               LWH/DataPointSet.h  \
  LWH/AIHistogramFactory.h        LWH/AIDataPointSetFactory.h     LWH/Histogram1D.h               LWH/DataPointSetFactory.h  \
  LWH/AITree.h                    LWH/AIMeasurement.h             LWH/HistogramFactory.h          LWH/Measurement.h  \
  LWH/AITreeFactory.h             LWH/AnalysisFactory.h           LWH/ManagedObject.h
endif
