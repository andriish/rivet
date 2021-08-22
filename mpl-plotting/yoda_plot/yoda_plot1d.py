import matplotlib.pyplot as plt
import numpy as np
import yoda
from collections.abc import Sequence
from .format_axis import format_axis, _get_axis_kw


def plot_hist(hists, plot_ref=True, ax=None, error_bars=True,
              colors=None, line_styles=['-', '--', '-.', ':'], legend=False, **kwargs):
    """Create a histogram plot. 

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
    if not isinstance(hists, Sequence):  # Convert single hist object to list
        hists = [hists]
    # Convert all histogram objects to Scatter2D
    hists = [h.mkScatter() if not isinstance(h, yoda.Scatter2D) else h for h in hists]
    if ax is None:
        ax = plt.gca()

    if isinstance(error_bars, bool):  # Convert to list of bool vals
        error_bars = [error_bars] * len(hists)
    if colors is None:  # Use default mpl prop cycle vals
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    line_properties = LineProperties(colors, line_styles)

    # Define useful variables
    x_points = (hists[0].xMins() + hists[0].xMaxs())/2
    x_bins = np.append(hists[0].xMins(), hists[0].xMax())
    ref_yvals = hists[0].yVals()

    # Plot the reference histogram data
    if plot_ref:
        ax.hlines(ref_yvals, hists[0].xMins(),
                  hists[0].xMaxs(), 'k')  # Plot reference data as horizontal lines
        ax.plot(x_points, ref_yvals, 'ko')  # Plot black dot in the middle of line
        ref_errminus = [err[0] for err in hists[0].yErrs()]
        ref_errplus = [err[1] for err in hists[0].yErrs()]
        if error_bars[0]:
            ax.vlines(x_points, (ref_yvals - ref_errminus),
                      (ref_yvals + ref_errplus), 'k')

    # Plot MC histogram data
    for i, hist in enumerate(hists):
        if plot_ref and i == 0:
            continue
        color, linestyle = next(line_properties)
        hist_yvals = np.insert(hist.yVals(), 0, hist.yVals()[0])
        ax.plot(x_bins, hist_yvals, color, linestyle=linestyle, drawstyle='steps-pre',
                solid_joinstyle='miter', zorder=5+i)
        if error_bars[i]:
            errminus = [err[0] for err in hist.yErrs()]
            errplus = [err[1] for err in hist.yErrs()]
            ax.vlines(x_points, (hist.yVals() - errminus),
                      (hist.yVals() + errplus), color, zorder=5+i, linestyle=linestyle)

    if legend:
        legend_labels = []
        for i, _ in enumerate(hists):
            if i == 0:
                legend_labels.append('Data')
            else:
                legend_labels.append('mc{}'.format(i))
        ax.legend(legend_labels)

    # Format x and y axes
    x_axis, y_axis, _ = _get_axis_kw(kwargs)
    if x_axis.get('lim') is None:
        x_axis['lim'] = [min([h.xMin() for h in hists]), max([h.xMax() for h in hists])]
    if y_axis.get('lim') is None and not y_axis.get('log'):
        y_axis['lim'] = 0
    format_axis('x', ax, **x_axis)
    format_axis('y', ax, **y_axis)

    return ax


def plot_ratio(hists, ax=None, error_bars=True, error_bands=False,
               colors=None, line_styles=['-', '--', '-.', ':'], legend=False, **kwargs):
    """Create a ratio plot of the reference and MC histograms.

    Parameters
    ----------
    hists : yoda histogram object or List[yoda histogram object]
        Yoda scatter objects that will be plotted where first element is reference data.
    ax : Axes object
        Axes object for plotting the histogram.
        If None, use `plt.gca()`.
    error_bars : bool or List[bool]
        Determines whether the error bars are shown on the plot. A single bool value corresponds to
        all the scatter objects. Otherwise, a list of bool values can be passed where the first element
        corresponds to the `ref_hist` object, and the remaining values correspond to the `mc_hists` list.
    error_bands : bool
        Determine whether the reference error is plotted as a band.
    colors : List[str]
        The list of colors to be used for plotting the scatter objects in the `mc_hists` list.
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
    if not isinstance(hists, Sequence):  # Convert single hist object to list
        hists = [hists]
    # Convert all histogram objects to Scatter2D
    hists = [h.mkScatter() if not isinstance(h, yoda.Scatter2D) else h for h in hists]
    if ax is None:
        ax = plt.gca()
    # Format x and y axes
    x_axis, y_axis, _ = _get_axis_kw(kwargs)
    if x_axis.get('lim') is None:
        x_axis['lim'] = [min([h.xMin() for h in hists]), max([h.xMax() for h in hists])]
    format_axis('x', ax, **x_axis)
    format_axis('y', ax, **y_axis)

    if isinstance(error_bars, bool):  # Convert to list of bool vals
        error_bars = [error_bars] * len(hists)
    if colors is None:  # Use default mpl prop cycle vals
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    line_properties = LineProperties(colors, line_styles)

    # Define useful variables
    x_points = (hists[0].xMins() + hists[0].xMaxs())/2
    x_bins = np.append(hists[0].xMins(), hists[0].xMax())
    ref_yvals = hists[0].yVals()

    # Plot the reference data
    ax.hlines(1, hists[0].xMin(), hists[0].xMax(), 'k', zorder=2)
    ax.plot(x_points, np.ones(len(x_points)), 'ko', zorder=3)
    if error_bars[0]:
        data_errminus = [err[0] for err in hists[0].yErrs()]
        data_errplus = [err[1] for err in hists[0].yErrs()]
        if error_bands:
            errbandminus = np.insert((ref_yvals - data_errminus)/ref_yvals, 0,
                                     ((ref_yvals - data_errminus)/ref_yvals)[0])
            errbandplus = np.insert((ref_yvals + data_errplus)/ref_yvals, 0,
                                    ((ref_yvals + data_errplus)/ref_yvals)[0])
            ax.fill_between(x_bins, errbandminus, errbandplus,
                            step='pre', alpha=0.5, zorder=0)
        else:
            ax.vlines(x_points, (ref_yvals - data_errminus)/ref_yvals,
                      (ref_yvals + data_errplus)/ref_yvals, 'k')

    # Plot the ratio data
    for i, hist in enumerate(hists[1:]):
        color, linestyle = next(line_properties)
        y_ratio = (np.insert(hist.yVals(), 0, hist.yVals()[0])
                   / np.insert(ref_yvals, 0, ref_yvals[0]))
        ax.plot(x_bins, y_ratio, color, linestyle=linestyle, drawstyle='steps-pre', zorder=1,
                solid_joinstyle='miter')
        if error_bars[i+1]:
            errminus = [err[0] for err in hist.yErrs()]
            errplus = [err[1] for err in hist.yErrs()]
            ax.vlines(x_points, (hist.yVals() - errminus)/ref_yvals,
                      (hist.yVals() + errplus)/ref_yvals, color, zorder=1, linestyle=linestyle)

    if legend:
        legend_labels = []
        for i, _ in enumerate(hists):
            if i == 0:
                legend_labels.append('Data')
            else:
                legend_labels.append('mc{}'.format(i))
        ax.legend(legend_labels)

    return ax


def plot_scatter1D(scatter1D, ax=None, colors=None, line_styles=['-', '--', '-.', ':'], legend=False, **kwargs):
    """Create a plot of Scatter1D objects. 

    Parameters
    ----------
    scatter1D : List[yoda.Scatter2D]
        Yoda 2Dscatter objects that will be plotted as MC data.
    ax : Axes object
        Axes object for plotting the histogram.
        If None, use `plt.gca()`.
    colors : List[str]
        The list of colors to be used for plotting the scatter objects in the `mc_hists` list.
        If None, the default matplotlib colors are used. 

    Returns
    -------
    ax
        The matplotlib axes object.
    """
    if not isinstance(scatter1D, Sequence):  # Convert single hist object to list
        scatter1D = [scatter1D]
    # Convert all histogram objects to Scatter2D
    scatter1D = [h.mkScatter() if not isinstance(h, yoda.Scatter1D) else h for h in scatter1D]
    if ax is None:
        ax = plt.gca()

    if colors is None:  # Use default mpl prop cycle vals
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    line_properties = LineProperties(colors, line_styles)

    # Plot the MC histogram data
    for i, mc in enumerate(scatter1D):
        # Cycle through colors for MC hists
        color, linestyle = next(line_properties)
        ax.hlines(mc.points()[0].val(1), 0, 1, color, linestyle=linestyle)
        ax.set_xlim(0, 1)
        ax.xaxis.set_visible(False)

    if legend:
        legend_labels = ['mc{}'.format(i+1) for i in range(len(scatter1D))]
        ax.legend(legend_labels)

    return ax


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
