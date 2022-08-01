import matplotlib.pyplot as plt
import numpy as np
import yoda
from collections.abc import Sequence
#from .format_axis import format_axis, _get_axis_kw

def getListsScatter1D(scatter1D, ax=None, colors=None, line_styles=['-', '--', '-.', ':'], legend=False, **kwargs):
  return """\nscatter1"""

def getListsHist(hists, plot_ref=True, ax=None, error_bars=True,
              colors=None, fnamedata=None, line_styles=['-', '--', '-.', ':'], legend=False, **kwargs):
    """
    Parameters
    ----------
    hists : yoda histogram object or List[yoda histogram object]
        Yoda scatter objects that will be plotted.
    plot_ref : Bool
        Determines whether the first index in `hists` is plotted as reference data.
    ax : Axes object
        Axes object for plotting the histogram.
        If None, use `plt.gca()`.
    error_bars : bool or List[bool]
        Determines whether the error bars are shown on the plot. A single bool value corresponds to
        all the scatter objects. Otherwise, a list of bool values can be passed where each element
        corresponds to the Scatter2D objects in `hists`.
    colors : List[str]
        The list of colors to be used for plotting the scatter objects in the `hists` list.
        If None, the default matplotlib colors are used. 
    line_styles : List[str]
        The list of line styles to be used in the `mc_hists` list.
    legend : Bool
        Determines whether to show a legend

    Other Parameters
    ----------------
    xlabel, ylabel : str, optional
        Axis label.
    xlim, ylim : tuple, optional
        Lower and upper limits of axis, in that order.
    logx, logy : bool, optional
        If True, set the scale of the axis to log scale. By default False
    xmajor_ticks, ymajor_ticks : int, optional
        Digit of the major ticks, by default None
    xminor_ticks, yminor_ticks : int, optional
        Number of minor ticks, by default None
    xcustom_major_ticks, ycustom_major_ticks : Sequence, optional
        Location of tick, followed by the corresponding tick label.
    xcustom_minor_ticks, ycustom_minor_ticks : Sequence[float], optional
        A list of the locations of minor ticks at arbitrary positions. 
    plot_xticklabels, plot_yticklabels : bool, optional
        If False, do not plot any tick labels, by default True
    **kwargs : Any
        Additional parameters passed to `ax.set`, e.g., `title`.

    Returns
    -------
    ax
        The matplotlib axes object.
    """

    #mplCommand = """""""
    listsWithData = {}
    legend_handles = []
    if not isinstance(hists, Sequence):  # Convert single hist object to list
        hists = [hists]
    # Convert all histogram objects to Scatter2D
    hists = [h.mkScatter() if not isinstance(h, yoda.Scatter2D) else h for h in hists]

    # Format x and y axes
    #x_axis, y_axis, _ = _get_axis_kw(kwargs)
    #if x_axis.get('lim') is None:
    #    x_axis['lim'] = [min([h.xMin() for h in hists]), max([h.xMax() for h in hists])]
    #if y_axis.get('lim') is None and not y_axis.get('log'):
    #    y_axis['lim'] = 0
    #format_axis('x', ax, **x_axis)
    #format_axis('y', ax, **y_axis)

    if isinstance(error_bars, bool):  # Convert to list of bool vals
        error_bars = [error_bars] * len(hists)
    if colors is None:  # Use default mpl prop cycle vals
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    line_properties = LineProperties(colors, line_styles)

    # Define useful variables
    x_points = (hists[0].xMins() + hists[0].xMaxs())/2
    x_bins = np.append(hists[0].xMins(), hists[0].xMax())
    ref_yvals = hists[0].yVals()
    ref_yerrs = hists[0].yErrs()
    lists = f"""
x_points  = {[val for val in x_points]}
x_bins    = {[val for val in x_bins]}
ref_yvals = {[val for val in ref_yvals]}
xMins     = {[val for val in hists[0].xMins()]}
xMaxs     = {[val for val in hists[0].xMaxs()]}
ref_errminus = {[err[0] for err in ref_yerrs]}
ref_errplus = {[err[1] for err in ref_yerrs]}
ref_errbarsdown = [ref_yvals[i] - ref_errminus[i] for i in range(len(ref_yvals))]
ref_errbarsup = [ref_yvals[i] + ref_errplus[i] for i in range(len(ref_yvals))]
"""

    # Plot the reference histogram data
    mplCommand = """# plot the data on axes\n"""
    fnamedata = fnamedata.split('py')[0]
    if plot_ref:
      mplCommand += f"""
# ref data
data_dots, = ax.plot(dataf.x_points, dataf.ref_yvals, 'ko', label='Data')  # Plot black dot in the middle of line
data_hlines = ax.hlines(dataf.ref_yvals, dataf.xMins, dataf.xMaxs, 'k')  # Plot reference data as horizontal lines
"""
      legend_handles += ['(data_dots, data_hlines)']
      if error_bars[0]: mplCommand += """ax.vlines(dataf.x_points, (dataf.ref_errbarsdown),(dataf.ref_errbarsup), 'k')""" 
    mplCommand += """\n
# histograms data from input yoda files"""
    # Plot MC histogram data
    for i, hist in enumerate(hists):
        if plot_ref and i == 0:
            continue
        color, linestyle = next(line_properties)
        hist_yvals = np.insert(hist.yVals(), 0, hist.yVals()[0])
        lists += f"""
mc{i}_yvals = {[val for val in hist_yvals]}
"""
        mplCommand += f"""
mc{i}, = ax.plot(dataf.x_bins, dataf.mc{i}_yvals, color='{color}', linestyle='{linestyle}', drawstyle='steps-pre',
        solid_joinstyle='miter', zorder={5+i})
"""
        legend_handles += [f'mc{i}']
        if error_bars[i]:
            errminus = [err[0] for err in hist.yErrs()]
            errplus = [err[1] for err in hist.yErrs()]
            mc_errbarsup = hist.yVals() - np.array(errminus)
            mc_errbarsdown = hist.yVals() + np.array(errplus)
            lists += f"""
mc{i}_errbarsup = {[val for val in mc_errbarsup]}
mc{i}_errbarsdown = {[val for val in mc_errbarsdown]}
"""
            mplCommand += f"""
ax.vlines(dataf.x_points, dataf.mc{i}_errbarsdown, dataf.mc{i}_errbarsup,
      color='{color}', zorder={5+i}, linestyle='{linestyle}')
"""
    mplCommand += f"""
handles = [{', '.join(legend_handles)}]
"""
# TODO ; set legend in main writePyScript function
#  if legend:
#      legend_labels = []
#      for i, _ in enumerate(hists):
#          if i == 0:
#              legend_labels.append('Data')
#          else:
#              legend_labels.append('mc{}'.format(i))
#      mplCommand += """\nax.legend(legend_labels)"""
    return mplCommand, lists


### 


def getListsRatio(hists, ax=None, error_bars=True, error_bands=False,
               colors=None, line_styles=['-', '--', '-.', ':'], legend=False, **kwargs):
  lists = """\n\n # lists for ratio plot"""
  mplCommand = """\n# Plot the ratio""" 
  if not isinstance(hists, Sequence):  # Convert single hist object to list
    hists = [hists]
  # Convert all histogram objects to Scatter2D
  hists = [h.mkScatter() if not isinstance(h, yoda.Scatter2D) else h for h in hists]

  if isinstance(error_bars, bool):  # Convert to list of bool vals
      error_bars = [error_bars] * len(hists)
  if colors is None:  # Use default mpl prop cycle vals
      colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
  line_properties = LineProperties(colors, line_styles)

  # Define useful variables
  #x_points = (hists[0].xMins() + hists[0].xMaxs())/2
  #x_bins = np.append(hists[0].xMins(), hists[0].xMax())
  ref_yvals = hists[0].yVals()

  # Plot the reference data
  mplCommand += f"""
{ax}.hlines(1, {hists[0].xMin()}, {hists[0].xMax()}, 'k', zorder=2)
{ax}.plot(dataf.x_points, [1.] * len(dataf.x_points), 'ko', zorder=3)
""" 
  if error_bars[0]:
    data_errminus = [err[0] for err in hists[0].yErrs()]
    data_errplus = [err[1] for err in hists[0].yErrs()]
    if error_bands:
      errbandminus = np.insert((ref_yvals - data_errminus)/ref_yvals, 0,
                               ((ref_yvals - data_errminus)/ref_yvals)[0])
      errbandplus = np.insert((ref_yvals + data_errplus)/ref_yvals, 0,
                              ((ref_yvals + data_errplus)/ref_yvals)[0])
      lists += f"""
ratio_errbandup = {[val for val in errbandplus]}
ratio_errbandminus = {[val for val in errbandminus]}
"""
      mplCommand += f"""
{ax}.fill_between(dataf.x_bins, dataf.ratio_errbandminus, dataf.ratio_errbandplus,
                step='pre', alpha=0.5, zorder=0)
"""
    else:
      lists += f"""
ratio_errbandup = {[val for val in (ref_yvals - data_errminus)/ref_yvals]}
ratio_errbandminus = {[val for val in (ref_yvals + data_errplus)/ref_yvals]} 
"""
      mplCommand += f"""
{ax}.vlines(dataf.x_points, dataf.ratio_errbandup, dataf.ratio_errbandminus, 'k')
"""  

  # Plot the ratio data
  for i, hist in enumerate(hists[1:]):
    color, linestyle = next(line_properties)
    y_ratio = (np.insert(hist.yVals(), 0, hist.yVals()[0])
            / np.insert(ref_yvals, 0, ref_yvals[0]))
    lists += f"""
ratio_y_mc{i} = {[val for val in y_ratio]}
"""
    mplCommand += f"""
{ax}.plot(dataf.x_bins, dataf.ratio_y_mc{i}, color='{color}', linestyle='{linestyle}', drawstyle='steps-pre', zorder=1,
         solid_joinstyle='miter')
"""
    print(error_bars[i+1])
    if error_bars[i+1]:
      errminus = [err[0] for err in hist.yErrs()]
      errplus = [err[1] for err in hist.yErrs()]
      lists += f"""
ratio_y_mc{i}_errsdown = {[val for val in (hist.yVals() - errminus)/ref_yvals]}
ratio_y_mc{i}_errsup   = {[val for val in (hist.yVals() + errplus)/ref_yvals]}
"""
      mplCommand += f"""
{ax}.vlines(dataf.x_points, dataf.ratio_y_mc{i}_errsdown, dataf.ratio_y_mc{i}_errsup,
                color='{color}', zorder=1, linestyle='{linestyle}')

"""
  # see mpl_tools/yoda_plot/yoda_plot1d.py -> plot_ratio function
   
  return mplCommand, lists 


####


class LineProperties:

    def __init__(self, colors, linestyles):
        self.colors = colors
        self.linestyles = linestyles
        self.color_index = 0
        self.style_index = 0

    def __iter__(self):
        return self

    def __next__(self):
        vals = (self.colors[self.color_index], self.linestyles[self.style_index])
        self.color_index += 1
        if self.color_index == len(self.colors):
            self.color_index = 0
            self.style_index += 1
            if self.style_index == len(self.linestyles):
                self.style_index = 0
        return vals
