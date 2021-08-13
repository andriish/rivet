import matplotlib.pyplot as plt
import numpy as np

#TODO: implement `format_axes` from yoda_plot2d.py when it's complete

def plot_hist(ref_hist, mc_hists=None, ax=None, ErrorBars=True, colors=None):
    """Create a histogram plot. 

    Parameters
    ----------
    ref_hist : yoda.Scatter2D
        Yoda scatter that will be plotted as reference data.
    mc_hists : List[yoda.Scatter2D]
        Yoda scatter objects that will be plotted as MC data.
    ax : Axes object
        Axes object for plotting the histogram.
        If None, use `plt.gca()`.
    ErrorBars : bool or List[bool]
        Determines whether the error bars are shown on the plot. A single bool value corresponds to
        all the scatter objects. Otherwise, a list of bool values can be passed where the first element
        corresponds to the `ref_hist` object, and the remaining values correspond to the `mc_hists` list.
    colors : List[str]
        The list of colors to be used for plotting the scatter objects in the `mc_hists` list.
        If None, the default matplotlib colors are used. 

    Returns
    -------
    ax
        The matplotlib axes object.
    """
    if ax is None:
        ax = plt.gca()
    if isinstance(ErrorBars, bool):  # Convert to list of bool vals
        if mc_hists is not None:
            ErrorBars = [ErrorBars] * (len(mc_hists)+1)
        else:
            ErrorBars = [ErrorBars]
    if colors is None:  # Use default mpl prop cycle vals
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    x_points = (ref_hist.xMins() + ref_hist.xMaxs())/2
    x_bins = np.append(ref_hist.xMins(), ref_hist.xMax())
    data_yVals = ref_hist.yVals()

    # Plot the reference histogram data
    ax.hlines(data_yVals, ref_hist.xMins(),
              ref_hist.xMaxs(), 'k')  # Plot reference data as horizontal lines
    ax.plot(x_points, data_yVals, 'ko')  # Plot black dot in the middle of line
    data_errminus = [err[0] for err in ref_hist.yErrs()]
    data_errplus = [err[1] for err in ref_hist.yErrs()]
    if ErrorBars[0]:
        ax.vlines(x_points, (data_yVals - data_errminus),
                  (data_yVals + data_errplus), 'k')

    # Plot the MC histogram data
    if mc_hists is not None:
        for i, mc in enumerate(mc_hists):
            # Cycle through colors for MC hists
            color = colors[i % len(colors)]
            y_mc = np.insert(mc.yVals(), 0, mc.yVals()[0])
            ax.plot(x_bins, y_mc, color, drawstyle='steps-pre',
                    solid_joinstyle='miter', zorder=5+i)  # Plot MC hist data
            if ErrorBars[i+1]:  # Plot MC error bars by default
                mc_errminus = [err[0] for err in mc.yErrs()]
                mc_errplus = [err[1] for err in mc.yErrs()]
                ax.vlines(x_points, (mc.yVals() - mc_errminus),
                          (mc.yVals() + mc_errplus), color, zorder=5+i)
    return ax


def plot_ratio(ref_hist, mc_hists=None, ax=None, ErrorBars=True, ErrorBands=False, colors=None):
    """Create a ratio plot of the reference and MC histograms.

    Parameters
    ----------
    ref_hist : yoda.Scatter2D
        Yoda scatter that will be plotted as reference data.
    mc_hists : List[yoda.Scatter2D]
        Yoda scatter objects that will be plotted as MC data.
    ax : Axes object
        Axes object for plotting the histogram.
        If None, use `plt.gca()`.
    ErrorBars : bool or List[bool]
        Determines whether the error bars are shown on the plot. A single bool value corresponds to
        all the scatter objects. Otherwise, a list of bool values can be passed where the first element
        corresponds to the `ref_hist` object, and the remaining values correspond to the `mc_hists` list.
    colors : List[str]
        The list of colors to be used for plotting the scatter objects in the `mc_hists` list.
        If None, the default matplotlib colors are used. 

    Returns
    -------
    ax
        The matplotlib axes object.
    """
    if ax is None:
        ax = plt.gca()
    if isinstance(ErrorBars, bool):  # Convert to list of bool vals
        if mc_hists is not None:
            ErrorBars = [ErrorBars] * (len(mc_hists)+1)
        else:
            ErrorBars = [ErrorBars]
    if colors is None:  # Use default mpl prop cycle vals
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    x_points = (ref_hist.xMins() + ref_hist.xMaxs())/2
    x_bins = np.append(ref_hist.xMins(), ref_hist.xMax())
    data_yVals = ref_hist.yVals()

    ax.hlines(1, ref_hist.xMin(), ref_hist.xMax(), 'k', zorder=2)
    ax.plot(x_points, np.ones(len(x_points)), 'ko', zorder=3)
    if ErrorBars[0]:
        data_errminus = [err[0] for err in ref_hist.yErrs()]
        data_errplus = [err[1] for err in ref_hist.yErrs()]
        if ErrorBands:
            errbandminus = np.insert((data_yVals - data_errminus)/data_yVals, 0,
                                     ((data_yVals - data_errminus)/data_yVals)[0])
            errbandplus = np.insert((data_yVals + data_errplus)/data_yVals, 0,
                                    ((data_yVals + data_errplus)/data_yVals)[0])
            ax.fill_between(x_bins, errbandminus, errbandplus,
                            step='pre', alpha=0.5, zorder=0)
        else:
            ax.vlines(x_points, (data_yVals - data_errminus)/data_yVals,
                      (data_yVals + data_errplus)/data_yVals, 'k')

    # Plot MCs histogram data
    if mc_hists is not None:
        for i, mc in enumerate(mc_hists):
            color = colors[i % len(colors)]
            y_ratio = (np.insert(mc.yVals(), 0, mc.yVals()[0])
                       / np.insert(data_yVals, 0, data_yVals[0]))
            ax.plot(x_bins, y_ratio, color, drawstyle='steps-pre', zorder=1,
                    solid_joinstyle='miter')
            if ErrorBars[i+1]:  # Plot MC error bars by default
                mc_errminus = [err[0] for err in mc.yErrs()]
                mc_errplus = [err[1] for err in mc.yErrs()]
                ax.vlines(x_points, (mc.yVals() - mc_errminus)/data_yVals,
                          (mc.yVals() + mc_errplus)/data_yVals, color, zorder=1)
    return ax

def plot_scatter1D(scatter1D, ax=None, colors=None):
    """Create a plot of scatter1D objects. 

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
    if ax is None:
        ax = plt.gca()
    if colors is None:  # Use default mpl prop cycle vals
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    x_bins = np.append(scatter1D[0].xMins(), scatter1D[0].xMax())

    # Plot the MC histogram data
    for i, mc in enumerate(scatter1D):
        # Cycle through colors for MC hists
        color = colors[i % len(colors)]
        ax.hlines(mc.points()[0].val(1), 0, 1, color)
        ax.set_xlim(0, 1)
        ax.xaxis.set_visible(False)
    return ax

