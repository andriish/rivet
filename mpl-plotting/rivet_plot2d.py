"""Create a rivet-style plot with 2D histograms."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yoda


def _plot_projection(hist_data, fig, hist_features, plot_features):
    axes = mpl.gridspec.gridspec.GridSpec(ncols=len(hist_data), nrows=2 if plot_features['RatioPlot'] else 1)
    colorbar_tick_format = mpl.ticker.ScalarFormatter(useMathText=True)
    colorbar_tick_format.set_powerlimits(
        (-plt.rcParams.get('axes.formatter.min_exponent'), plt.rcParams.get('axes.formatter.min_exponent'))
    )
    
    for h, ax in zip(hist_data, axes):
        x_edges, y_edges = np.meshgrid(h.xEdges(), h.yEdges())
        im = ax.pcolormesh(x_edges, y_edges, h.zVals(asgrid=True).T, norm=mpl.colors.Normalize(vmin=h.zMin(), vmax=h.zMax()))
        fig.colorbar(im, ax=ax, format=colorbar_tick_format)
        ax.set_title(h['Title'])

    return axes


def _plot_surface(hist_data, fig, hist_features, plot_features):
    axes = mpl.gridspec.gridspec.GridSpec(ncols=len(hist_data), nrows=2 if plot_features['RatioPlot'] else 1)
    
    return axes


def _plot_contour(hist_data, fig, hist_features, plot_features):
    # Set nrows=1 here when implementing
    raise NotImplementedError('2DType "contour" has not been implemented yet. '
        'Use "projection" or "surface" instead.')


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
    # TODO these 2 lines will be removed when axes are always created inside the plotting function
    fig = axes.get_figure()
    fig.clf()
    if plot_features.get('2DType', 'projection') == 'projection':
        _plot_projection(hist_data, fig, hist_features, plot_features)
    elif plot_features.get('2DType') == 'surface':
        _plot_surface(hist_data, fig, hist_features, plot_features)
    elif plot_features.get('2DType') == 'contour':
        _plot_contour(hist_data, fig, hist_features, plot_features)
    else:
        raise NotImplementedError('Expected 2DType to be "projection", "surface" or "contour" but got {}' % plot_features.get('2DType'))
