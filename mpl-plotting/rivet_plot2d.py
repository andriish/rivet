"""Create a rivet-style plot with 2D histograms."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yoda


def _plot_projection(hist_data, fig, hist_features, plot_features):
    gs = fig.add_gridspec(ncols=len(hist_data), nrows=2 if plot_features.get('RatioPlot', True) else 1)
    colorbar_tick_format = mpl.ticker.ScalarFormatter(useMathText=True)
    colorbar_tick_format.set_powerlimits(
        (-plt.rcParams.get('axes.formatter.min_exponent'), plt.rcParams.get('axes.formatter.min_exponent'))
    )
    
    for i, (h, settings) in enumerate(zip(hist_data, hist_features)):
        ax = fig.add_subplot(gs[0, i])
        x_edges, y_edges = np.meshgrid(h.xEdges(), h.yEdges())
        # TODO: should colorbar have min and max for all of them instead? In this case, only one colorbar is used.
        im = ax.pcolormesh(x_edges, y_edges, h.zVals(asgrid=True).T, norm=mpl.colors.Normalize(vmin=h.zMin(), vmax=h.zMax()))
        cbar = fig.colorbar(im, ax=ax, format=colorbar_tick_format)
        cbar.set_label(plot_features.get('ZLabel'))
        ax.set(xlabel=plot_features.get('XLabel'), ylabel=plot_features.get('YLabel'), title=settings['Title'])

    if plot_features.get('RatioPlot', True):
        ref_h, ref_settings = hist_data[0], hist_features[0]
        for i, h in enumerate(hist_data[1:]):
            ax = fig.add_subplot(gs[1, i+1])
            x_edges, y_edges = np.meshgrid(h.xEdges(), h.yEdges())
            im = ax.pcolormesh(x_edges, y_edges, (h.zVals(asgrid=True) / ref_h.zVals(asgrid=True)).T)
            cbar = fig.colorbar(im, ax=ax, format=colorbar_tick_format)
            cbar.set_label(plot_features.get('ZLabel'))
            ax.set(xlabel=plot_features.get('XLabel'), ylabel=plot_features.get('YLabel'), 
                title=settings.get('RatioPlotTitle', '{}/{}'.format(settings['Title'], ref_settings['Title'])))
    fig.tight_layout()


def _plot_surface(hist_data, fig, hist_features, plot_features):
    axes = mpl.gridspec.gridspec.GridSpec(ncols=len(hist_data), nrows=2 if plot_features.get('RatioPlot', True) else 1)
    
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
        TODO: currently only supports the "Title" setting but more should be added. 
    plot_features : dict
        Plot settings for the entire axes.

    Returns
    -------
    None
    """
    # TODO these 2 lines will be removed when axes are always created inside the plotting function
    fig = axes[0].get_figure()
    fig.clf()
    if plot_features.get('2DType', 'projection') == 'projection':
        _plot_projection(hist_data, fig, hist_features, plot_features)
    elif plot_features.get('2DType') == 'surface':
        _plot_surface(hist_data, fig, hist_features, plot_features)
    elif plot_features.get('2DType') == 'contour':
        _plot_contour(hist_data, fig, hist_features, plot_features)
    else:
        raise NotImplementedError('Expected 2DType to be "projection", "surface" or "contour" but got {}'.format(plot_features.get('2DType')))
