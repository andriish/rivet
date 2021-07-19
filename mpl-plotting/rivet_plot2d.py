"""Module creates a rivet-style plot with 2D histograms."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yoda


def _plot_2Dhist(hist_data, axes, hist_features, plot_features):
    """Plot 2D histogram in hist_data based 

    Parameters
    ----------
    hist_data : list[yoda.Histo2D | yoda.Profile2D | yoda.Scatter3D]
        All histograms that will be plotted.
    axes : matplotlib.axes.Axes
        The axes object in which the histogram(s) will be plotted.
        TODO axes will have to be created inside this function later.
    hist_features : dict
        Plot settings for each histogram.
    plot_features : dict
        Plot settings for the entire axes.

    Returns 
    -------
    None
    """
    fig = axes.get_figure()
    colorbar_tick_format = mpl.ticker.ScalarFormatter(useMathText=True)
    colorbar_tick_format.set_powerlimits((-plt.rcParams.get('axes.formatter.min_exponent'), plt.rcParams.get('axes.formatter.min_exponent')))

    for h, ax in zip(hist_data, axes):
        x_edges, y_edges = np.meshgrid(h.xEdges(), h.yEdges())
        im = ax.pcolormesh(x_edges, y_edges, h.zVals(asgrid=True).T, cmap='jet', norm=mpl.colors.Normalize(vmin=h.zMin(), vmax=h.zMax()))
        fig.colorbar(im, ax=ax, format=colorbar_tick_format)
        
        fig.colorbar(im, ax=ax, format=colorbar_tick_format)
        ax.set_title(h.title() if h.title() != '~' else h.name())
