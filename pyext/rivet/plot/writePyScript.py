import os

import yoda
import rivet

from rivet.plot.mpl_tools.yamlio import read_yamlfile
from rivet.plot.mpl_tools.mathtext_preprocessor import preprocess
from rivet.plot import writePyScript_fetchData

import matplotlib as mpl
import matplotlib.pyplot as plt

import math
import numpy
import re 

def mkCanvas(ratio, yoda_type):
  """ Return Python command to plot a basic canvas in matplotlib
  
  """

  # Create fig and axes
  if ratio and yoda_type == 'hist':
    return """\n
# create figure and axis objects
fig, (ax, ax_ratio) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': (2, 1)})"""
  else:
    return """\n
# create figure and axis objects
fig, ax = plt.subplots(1, 1)"""
  
def writePyScript1D(hist_data, hist_features, yaml_dicts, outdir, plot_name):
  """ Return Python commands to make a Rivet plot from the input parameters
  """

  mplCommand1D = ""
  
  # set style 
  mpl_stylename = yaml_dicts.get('style', 'rivet') + '.mplstyle'
  plot_style = rivet.findAnalysisPlotFile(os.path.join("plot", mpl_stylename)) 
  if not os.path.isfile(plot_style):
      raise NotImplementedError('Plot style file not found.')
  
  # place copy of plot style in out directory
  if not os.path.isfile(os.path.join(outdir, mpl_stylename)): os.system(f"cp {plot_style} {os.path.join(outdir, mpl_stylename)}")
 
  # set plot style 
  plt.style.use(plot_style)
  mplCommand1D += f"""\n#plot style \nplt.style.use('{os.path.join(outdir, mpl_stylename)}')"""
  
  plot_features = yaml_dicts.get('plot features', {})

  mplCommand1D += f"""\n
# plot metadata
ax_xLabel = r'{plot_features.get('XLabel')}'
ax_yLabel = r'{plot_features.get('YLabel')}'
ax_title  = r'{plot_features.get('Title')}'"""

  xscale = 'log' if (plot_features.get('LogX') ) else 'linear'
  yscale = 'log' if (plot_features.get('LogY', 1) ) else 'linear'
  mplCommand1D += f"""
ax_xScale = '{xscale}'
ax_yScale = '{yscale}'
"""
  yoda_type = 'hist' if isinstance(hist_data[0], yoda.Scatter2D) else 'scatter'

  # TODO there should be a cleaner way, but \textrm is generally
  # not working in math mode. Naively check if there is a regex match
  # to a \textrm command in math mode to catch some exceptions... 
  if 'Title' in plot_features.keys() and re.search("^.*\$.*textrm.*\$.*$", plot_features['Title']):
    plot_features['Title'] = plot_features['Title'].replace("\\textrm", "\\mathrm")
  if 'XLabel' in plot_features.keys() and re.search("^.*\$.*textrm.*\$.*$", plot_features['XLabel']):
    plot_features['XLabel'] = plot_features['XLabel'].replace("\\textrm", "\\mathrm")
  if 'YLabel' in plot_features.keys() and re.search("^.*\$.*textrm.*\$.*$", plot_features['YLabel']):
    plot_features['YLabel'] = plot_features['YLabel'].replace("\\textrm", "\\mathrm")

  ax_format = {}  # Stores the items for formatting the axes in a dict
  
  # Set plot lims
  if yoda_type == 'hist':
      XMin = plot_features.get('XMin', min([h.xMin() for h in hist_data]))
      XMax = plot_features.get('XMax', max([h.xMax() for h in hist_data]))
      ax_format['xlim'] = (XMin, XMax)

  # set maximum Y value from all hist datasets 
  if yoda_type == 'scatter':
      max_ymax = max([h.points()[0].val(1) for h in hist_data])
  else:
      max_ymax = max([max(h.yVals()) for h in hist_data])
  if plot_features.get('YMax') is not None:
      YMax = plot_features.get('YMax')
  elif plot_features.get('LogY', 1):
      # round off highest number in the histograms to next power of 10
      YMax = 10**(math.ceil(math.log10(max_ymax)))
  else:
      YMax = 1.1*max_ymax
  
  # Use minimum y value from all hist datasets
  if yoda_type == 'scatter':
      min_ymin = min([h.points()[0].val(1) for h in hist_data]) # TO DO -- where does this come from??
  else:
      min_ymin = min([min(h.yVals()) for h in hist_data])
  if plot_features.get('YMin') is not None:
      YMin = plot_features.get('YMin')
  elif plot_features.get('LogY', 1):
      # round off lowest number in the histograms to lower power of 10
      YMin = 10**(math.floor(math.log10(min_ymin))) if min_ymin !=0 else 2e-7*YMax # TODO: come up with a better solution to deal with min_ymin=0
  elif plot_features.get('ShowZero', 1):  # default ShowZero is True
      YMin = 0 if min_ymin > -1e-4 else 1.1*min_ymin
  else:
      YMin = (1.1*min_ymin if min_ymin < -1e-4 else 0 if min_ymin < 1e-4
              else 0.9*min_ymin)

  ax_format['ylim'] = (YMin, YMax)
  ax_format['logx'] = plot_features.get('LogX')
  ax_format['logy'] = plot_features.get('LogY', 1)

  if ax_format['xlim']: mplCommand1D += f"""\nxLims = {ax_format['xlim']}"""
  if ax_format['ylim']: mplCommand1D += f"""\nyLims = {ax_format['ylim']}"""
 
  # labels for the legend
  labels = []
  if plot_features.get('Legend', 1) and yoda_type == 'hist':
      labels = [hist_features[0].get('Title', 'Data')]
      for i, _ in enumerate(hist_data[1:]):
          labels.append(hist_features[i+1].get('Title', 'mc{}'.format(i+1)).replace(".yoda",""))
  if plot_features.get('Legend', 1) and yoda_type == 'scatter':
      labels = []
      for i, _ in enumerate(hist_data):
          labels.append(hist_features[i].get('Title', 'mc{}'.format(i+1)).replace(".yoda",""))
  if len(labels) > 0:
    mplCommand1D += f"""\n
# labels for in the legend
labels = {labels}"""

  mkRatio = plot_features.get('RatioPlot', 1) 
  mplCommand1D += mkCanvas(mkRatio, yoda_type)
  if mkRatio and yoda_type == 'hist': 
    # Ratio plot has y range of 0.5 to 1.5
    RatioPlotYMin = plot_features.get('RatioPlotYMin', 0.5)
    RatioPlotYMax = plot_features.get('RatioPlotYMax', 1.4999)  # Don't plot 1.5
    mplCommand1D += f"""
ax_ratio.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
ax_ratio.set_ylabel('{plot_features.get('RatioPlotYLabel', 'MC/Data')}')
ax_ratio.set_ylim({RatioPlotYMin}, {RatioPlotYMax})
""" 
    
  if ax_format['logy']:
    mplCommand1D += f"""
ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(
    base=10.0, subs=[i for i in np.arange(0, 1, 0.1)], numticks=np.inf))
"""
  ax_format['xmajor_ticks'] = plot_features.get('XMajorTickMarks')
  ax_format['ymajor_ticks'] = plot_features.get('YMajorTickMarks')
  ax_format['xminor_ticks'] = plot_features.get('XMinorTickMarks')
  ax_format['yminor_ticks'] = plot_features.get('YMinorTickMarks')
  ax_format['xcustom_major_ticks'] = plot_features.get('XMajorTickMarks')
  #if plot_features.get('nRatioTicks'): print("NRATIOTICKS: ", plot_features.get('nRatioTicks'))

  if plot_features.get('XCustomMajorTicks') is not None:
      ax_format['xcustom_major_ticks'] = plot_features.get('XMajorTickMarks')
      if plot_features.get('RatioPlot', 1):
        mplCommand1D += """\nax_ratio.set_xticks([], minor=True)"""

  ax_format['ycustom_major_ticks'] = plot_features.get('YMajorTickMarks')
  ax_format['ycustom_minor_ticks'] = plot_features.get('YMajorTickMarks')
  ax_format['ycustom_minor_ticks'] = plot_features.get('YMajorTickMarks')
  ax_format['plot_xticklabels'] = plot_features.get('PlotXTickLabels') 
  
  plot_errorbars = [h.get('ErrorBars', 1) for h in hist_features]
  colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
  mplCommand1D += f"""colors = {colors}"""

  dataOutPyName = plot_name.strip('/').replace('-', '_') + '__data.py'
  dataOutFile = open(os.path.join(outdir, dataOutPyName), "a")
  print(os.path.join(outdir, dataOutPyName))
  mplCommand1D += f"""\n
# the numerical data is stored in a separate file
sys.path.insert(1, './rivet-plots')
import {dataOutPyName.split('.py')[0].replace('-', '_').replace('/', '.')} as dataf""" 
  if yoda_type == 'scatter':
    mplCommand1D += getListsScatter1D(hist_data, "ax", colors, **ax_format) #TODO
  else:
    fetchHistData = rivet.plot.writePyScript_fetchData.getListsHist(hist_data, True, "ax", plot_errorbars, colors, dataOutPyName, **ax_format)
    mplCommand1D += fetchHistData[0]
    dataOutFile.write(fetchHistData[1])
    if plot_features.get('RatioPlot', 1):
      fetchRatioData = rivet.plot.writePyScript_fetchData.getListsRatio(hist_data, "ax_ratio", plot_errorbars, plot_features.get('ErrorBands'), colors)
      mplCommand1D += fetchRatioData[0]
      dataOutFile.write(fetchRatioData[1])
      
  dataOutFile.close()

  

  if plot_features.get('Legend', 1):
      if plot_features.get('LegendAlign') is None or plot_features.get('LegendAlign') == 'l':
          legend_pos = (plot_features.get('LegendXPos', 0.5),
                        plot_features.get('LegendYPos', 0.97))
          mplCommand1D += f"""
legend_pos = {legend_pos}
ax.legend(legend_handles, labels, loc='upper left', bbox_to_anchor=legend_pos)"""
      if plot_features.get('LegendAlign') == 'r':
          legend_pos = (plot_features.get('LegendXPos', 0.97),
                        plot_features.get('LegendYPos', 0.97))
          mplCommand1D += f"""
legend_pos = {legend_pos}

ax.legend(legend_handles, labels, loc='upper right', bbox_to_anchor=legend_pos,markerfirst=False)
"""

  # Set text labels on axes
  mplCommand1D += """\n# set plot metadata as defined above"""
  if plot_features.get('RatioPlot', 1) and yoda_type == 'hist':
      mplCommand1D += f"""\n
ax_ratio.set_xlabel(ax_xLabel)"""
  else:
      mplCommand1D += f"""\n
ax.set_xlabel(ax_xLabel)"""
  mplCommand1D += f"""
ax.set_ylabel(ax_yLabel, loc='top')
ax.set_title(ax_title, loc='left')
ax.set_xscale(ax_xScale)
ax.set_yscale(ax_yScale)"""


  # toggle x/y lims
  if ax_format['xlim']: mplCommand1D += """\nax.set_xlim(xLims)"""
  if ax_format['ylim']: mplCommand1D += """\nax.set_ylim(yLims)"""
  if plot_features.get('RatioPlot', 1) and yoda_type == 'hist': mplCommand1D += """\n\nfig.align_ylabels((ax, ax_ratio))"""

  mplCommand1D += f"""

# tick formatting
plt.rcParams['xtick.top'] = {plot_features.get('XTwosidedTicks', True)}
plt.rcParams['ytick.right'] = {plot_features.get('YTwosidedTicks', True)}
"""
  
  # temporary for debugging
  #mplCommand1D += """\nplt.show()"""
  mplCommand1D += f"""\nplt.savefig('{os.path.join(outdir, plot_name.strip('/')) + '.pdf'}')"""
  return mplCommand1D

def writePyScript2D():
  #TODO 
  pass

def writePyScript(yaml_file, plot_name, outdir):
  """ write an executable python script that only relies on YODA
      functionalities
  
  Parameters
  ----------
  yaml_file : TODO
  plot_name : TODO
  outdir    : TODO
  """
    
  outPyName = os.path.join(outdir, plot_name.strip('/')) + '.py'

  # contents to write out to .py script
  mplScript = """#! /usr/bin/env python"""

  mplScript +='''
"""
This Python script creates your Rivet plots and is implicitly called
when running rivet-mkhtml. It is written out explicitly to give the 
option to change the plot features, such as axis-labels, line colours 
or the position of the legend.
"""
'''

  # add imports
  mplScript += """
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import numpy as np
"""
  
  # Parse yaml file for histogram data
  if isinstance(yaml_file, str):  # If the file name of the yaml file is passed
      yaml_dicts = read_yamlfile(yaml_file)
      hist_data = _parse_yoda_hist(yaml_dicts)
  else:  # If the dictionary object is passed
      yaml_dicts = yaml_file
      hist_data = [
          hist_dict['yoda'].mkScatter()
          if not isinstance(hist_dict, (yoda.Scatter1D, yoda.Scatter2D, yoda.Scatter3D)) else hist_dict['yoda']
          for hist_dict in yaml_dicts['histograms'].values()
      ]
  
  hist_features = [val for val in yaml_dicts['histograms'].values()]
  output_filename = os.path.join(outdir, plot_name.strip('/'))
  _preprocess_text(yaml_dicts)

  # Ensure reference histogram is first in list since dicts are not ordered in Python 2, not needed for Python 3 
  # TODO remove since not needed for Py3?
  #for i, h in enumerate(hist_features):
  #    if h.get('IsRef'):
  #        if i != 0:  # Not necessary since the reference histogram is already the first element
  #            hist_data.insert(0, hist_data.pop(i))
  #            hist_features.insert(0, hist_features.pop(i))
  #        break

  if all(isinstance(h, yoda.Scatter1D) for h in hist_data) or all(isinstance(h, yoda.Scatter2D) for h in hist_data):
      mplScript += writePyScript1D(hist_data, hist_features, yaml_dicts, outdir, plot_name)

  elif all(isinstance(h, yoda.Scatter3D) for h in hist_data):
      plot_2Dhist(hist_data, hist_features, yaml_dicts, output_filename)
  else:
      print('Error with Class types:', [type(h) for h in hist_data])
      raise NotImplementedError('Class type cannot be plotted yet')
  
  with open(outPyName, "a") as f:
    f.write(mplScript)
    #f.write("import matplotlib.pyplot as plt")

  return(outPyName)
   
# ultimately, have something produced of this form:
  """
    #!/usr/bin/python
    import matplotlib.pyplot as plt
    from YODA import mkplots, mkCanvas

    style = 'rivet'
    if <more than two input files> or <found Ref data and at least one input .yoda>:
      ratio='True'
    else:
      ratio='False'
    fig,ax = yoda.mkCanvas(style, ratio='True') # --> this makes the canvas, tick size, etc. etc.
    
    ax.set_ylim(5, 20)
    ax.set_yscale('log')
    ax.set_xscale('lin')
    plt.rcParams['legend.loc'] = "upper right"

    # these lists need to contain explicit numerical values, otherwise can't
    # run it separately from rivet
    plotData = {'refData' : <list>, 'yoda1' : <list>, 'yoda2' : <list>}
    bins = <list>
    yoda.mkPlots(fig, ax, bins, plotData['refData'], color='black', type='Data')
    yoda.mkPlots(fig, ax, bins, plotData['yoda1'],    color='green', type='MC') 
    yoda.mkPlots(fig, ax, bins, plotData['yoda2'],    color='green', type='MC') 
    
    VECTORFORMAT = 'pdf' # command line option with rivet-mkhtml
    plt.savefig(f'/path/to/file.{VECTORFORMAT}')
  """

  
  pass  

def _parse_yoda_hist(yaml_dicts):
    """Read yoda string and return yoda object."""
    hist_data = []
    for hist_dict in yaml_dicts['histograms'].values():
        with io.StringIO(hist_dict['yoda']) as file_like:
            hist_data.append(yoda.readYODA(
                file_like, asdict=False)[0].mkScatter())
    return hist_data


def _preprocess_text(yaml_dicts):
    """Preprocess text to convert convenient HEP units and other symbols to mathtext."""
    if 'plot features' not in yaml_dicts:
        yaml_dicts['plot features'] = {}
    plot_features = yaml_dicts['plot features']
    for plot_property in ('Title', 'XLabel', 'YLabel', 'ZLabel', 'RatioPlotYLabel'):
        if plot_property in plot_features:
            plot_features[plot_property] = preprocess(plot_features[plot_property])
    for custom_major_ticks in ('XCustomMajorTicks', 'YCustomMajorTicks', 'ZCustomMajorTicks'):
        if custom_major_ticks in plot_features:
            plot_features[custom_major_ticks] = [preprocess(tick) for tick in plot_features[custom_major_ticks]]

    for hist_setting in yaml_dicts['histograms'].values():
        if 'Title' in hist_setting:
            hist_setting['Title'] = preprocess(hist_setting['Title'])

class AnyObject(object):
    """Necessary for custom legend handler."""
    pass


class AnyObjectHandler(object):
    """Creates custom legend handler for Data."""

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        patch = mpl.patches.Circle(
            (x0+width/2, y0+height/2), (2.5)**0.5, facecolor='black')
        handlebox.add_artist(patch)
        patch = mpl.patches.Rectangle(
            (-0.4, 3), width+0.8, 0.8, facecolor='black')
        handlebox.add_artist(patch)
        patch = mpl.patches.Rectangle(
            (width/2-0.4, 0), 0.8, height, facecolor='black')
        handlebox.add_artist(patch)
        return patch
