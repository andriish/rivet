"""Create a rivet-style plot with 2D histograms."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
# TODO
#   Only have 1 colorbar per row (?)
#   Implement ShowZero
#   Implment contour plot
#   Should the ratio plots touch the 2D histogram?
#   Spacing between plots and colorbar
def _configure_plot(ax, plot_features, hist_data):
    # TODO: optimize using e.g. sharex.
    # TODO: Different settings for ratio and non-ratio?
    ax.set_xlabel(plot_features.get('XLabel'))
    ax.set_ylabel(plot_features.get('YLabel'), loc='top')
    ax.set_title(plot_features.get('Title'), loc='left')

    # TODO: check that this is ok. Maybe more customization is needed.
    for axis_name in 'XY':
        min_lim = plot_features.get(axis_name + 'Min', min(getattr(h, axis_name.lower() + 'Min')() for h in hist_data))
        max_lim = plot_features.get(axis_name + 'Max', max(getattr(h, axis_name.lower() + 'Max')() for h in hist_data))
        ax.set(**{axis_name.lower() + 'lim': (min_lim, max_lim)})
    # TODO: zlim. Different when in 3D and when using color

    # Set log scale
    if plot_features.get('LogX'):
        ax.set_xscale('log')
        ax.xaxis.set_major_locator(mpl.ticker.LogLocator(numticks=np.inf))
    if plot_features.get('LogY'):
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(mpl.ticker.LogLocator(numticks=np.inf))
        ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(
            base=10.0, subs=[i for i in np.arange(0, 1, 0.1)], numticks=np.inf))
    # TODO: LogZ

    # Set tick marks frequency
    if plot_features.get('XMajorTickMarks') is not None and not plot_features.get('LogX'):
        base = plot_features.get('XMajorTickMarks')*10**(int(np.log10(XMax))-1)
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base))
    if plot_features.get('YMajorTickMarks') is not None and not plot_features.get('LogY'):
        base = plot_features.get('YMajorTickMarks')*10**(int(np.log10(YMax))-1)
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base))

    if plot_features.get('XMinorTickMarks') is not None and not plot_features.get('LogX'):
        ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(
            1+plot_features.get('XMinorTickMarks')))
    if plot_features.get('YMinorTickMarks') is not None and not plot_features.get('LogY'):
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(
            1+plot_features.get('YMinorTickMarks')))
    # TODO: Z tick marks

    # Add custom ticks
    if plot_features.get('XCustomMajorTicks') is not None:
        ax.set_xticks(plot_features.get('XCustomMajorTicks')[::2])
        ax.set_xticklabels(plot_features.get('XCustomMajorTicks')[1::2])
        ax.set_xticks([], minor=True)  # Turn off minor xticks
    if plot_features.get('YCustomMajorTicks') is not None:
        ax.set_yticks(plot_features.get('YCustomMajorTicks')[::2])
        ax.set_yticklabels(plot_features.get('YCustomMajorTicks')[1::2])
        ax.set_yticks([], minor=True)  # Turn off minor yticks
    if plot_features.get('XCustomMinorTicks') is not None:
        ax.set_xticks(plot_features.get('XCustomMinorTicks'), minor=True)
    if plot_features.get('YCustomMinorTicks') is not None:
        ax.set_yticks(plot_features.get('YCustomMinorTicks'), minor=True)
    if plot_features.get('PlotXTickLabels') == 0:
        ax.set_xticklabels([])
    # TODO: Z tick labels


def _plot_projection(hist_data, fig, hist_features, plot_features):
    plot_ratio = plot_features.get('RatioPlot', True) and len(hist_data) > 1
    ncols = len(hist_data) + 1
    nrows = 2 if plot_ratio else 1
    # TODO Try different values for size of colorbar
    width_ratios = np.array([1] * len(hist_data) + [0.0526])
    fig.set_size_inches(np.sum(width_ratios) * plt.rcParams['figure.figsize'][0], plt.rcParams['figure.figsize'][1] * ncols)
    
    gs = fig.add_gridspec(ncols=ncols, nrows=nrows, width_ratios=width_ratios)
    colorbar_tick_format = mpl.ticker.ScalarFormatter(useMathText=True)
    colorbar_tick_format.set_powerlimits(
        (-plt.rcParams.get('axes.formatter.min_exponent'), plt.rcParams.get('axes.formatter.min_exponent'))
    )
    
    zmin, zmax = min(h.zMin() for h in hist_data), max(h.zMax() for h in hist_data)

    for i, (h, settings) in enumerate(zip(hist_data, hist_features)):
        ax = fig.add_subplot(gs[0, i])
        x_edges, y_edges = np.meshgrid(h.xEdges(), h.yEdges())
        # TODO: should colorbar have min and max for all of them instead? In this case, only one colorbar is used.
        im = ax.pcolormesh(x_edges, y_edges, h.zVals(asgrid=True).T, norm=mpl.colors.Normalize(vmin=zmin, vmax=zmax))
        _configure_plot(ax, plot_features, hist_data)
    cbar = fig.colorbar(im, cax=fig.add_subplot(gs[0, -1]), format=colorbar_tick_format)
    cbar.set_label(plot_features.get('ZLabel'))

    if plot_features.get('RatioPlot', True):
        ref_h, ref_settings = hist_data[0], hist_features[0]
        zmin = zmax = np.nan
        for i, h in enumerate(hist_data[1:]):
            ax = fig.add_subplot(gs[1, i+1])
            x_edges, y_edges = np.meshgrid(h.xEdges(), h.yEdges())
            
            z_ratio = (h.zVals(asgrid=True) / ref_h.zVals(asgrid=True)).T
            z_ratio[~np.isfinite(z_ratio)] = np.nan
            zmin_new, zmax_new = np.nanmin(z_ratio), np.nanmax(z_ratio)
            zmin, zmax = np.nanmin([zmin, zmin_new]), np.nanmax([zmax, zmax_new])            
            im = ax.pcolormesh(x_edges, y_edges, z_ratio)
            _configure_plot(ax, plot_features, hist_data)

        cbar = fig.colorbar(im, cax=fig.add_subplot(gs[1, -1]), format=colorbar_tick_format)
        cbar.set_label(plot_features.get('ZLabel'))
            

def _plot_surface(hist_data, fig, hist_features, plot_features):
    """Plot 3D plots with a surface representing histograms.

    Parameters
    ----------
    hist_data : list[yoda.Histo2D | yoda.Profile2D | yoda.Scatter3D]
        All the YODA histograms that will be plotted. Histograms will be plotted in separate axes.
    fig : plt.figure.Figure
        Matplotlib figure object in which all the subplots will be plotted in.
    hist_features : list[dict]
        Plot settings for each histogram, in the same order as hist_data
    plot_features : dict
        Settings that will be applied to the entire figure or each axes.
    """
    plot_ratio = plot_features.get('RatioPlot', True) and len(hist_data) > 1
    ncols = len(hist_data)
    nrows = 2 if plot_ratio else 1
    gs = fig.add_gridspec(ncols=ncols, nrows=nrows)

    for i, (h, settings) in enumerate(zip(hist_data, hist_features)):
        ax = fig.add_subplot(gs[0, i], projection='3d')
        ax.plot_surface(h.xVals(asgrid=True), h.yVals(asgrid=True), h.zVals(asgrid=True), 
            norm=mpl.colors.Normalize(vmin=h.zMin(), vmax=h.zMax()), cmap=plt.rcParams['image.cmap'])
        # TODO Move these kinds of settings to separate function, similar to rivet_plot._create_plot 
        ax.set(xlabel=plot_features.get('XLabel'), ylabel=plot_features.get('YLabel'), title=settings['Title'])

    if plot_features.get('RatioPlot', True):
        ref_h, ref_settings = hist_data[0], hist_features[0]
        for i, h in enumerate(hist_data[1:]):
            ax = fig.add_subplot(gs[1, i+1], projection='3d')
            z_ratio = h.zVals(asgrid=True) / ref_h.zVals(asgrid=True)
            z_ratio[~np.isfinite(z_ratio)] = np.nan
            ax.plot_surface(h.xVals(asgrid=True), h.yVals(asgrid=True), z_ratio,
                cmap=plt.rcParams['image.cmap'], vmin=np.min(z_ratio[np.isfinite(z_ratio)]), vmax=np.max(z_ratio[np.isfinite(z_ratio)]))
            # TODO Move these kinds of settings to separate function, similar to rivet_plot._create_plot 
            ax.set(xlabel=plot_features.get('XLabel'), ylabel=plot_features.get('YLabel'), 
                title=settings.get('RatioPlotTitle', '{}/{}'.format(settings['Title'], ref_settings['Title'])))


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
    fig : Figure
        Matplotlib Figure object containing all plots.
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
    
    fig.suptitle(plot_features.get('Title'))
    # _configure_plot(fig)

    return fig
