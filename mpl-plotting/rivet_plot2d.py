"""Create a rivet-style plot with 2D histograms."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
# TODO
#   Only have 1 colorbar per row (?)
#     Spacing between plots and colorbar
#   Implement ShowZero
#   Implement contour plot
#   Should the ratio plots touch the 2D histogram?
#   Make separate mplstyle for 2d plots?

def _configure_plot(ax, plot_features, hist_data, hist_settings=None):
    # TODO: Different settings for ratio and non-ratio?
    if hist_settings is None:
        hist_settings = {}
    ax.set_title(hist_settings.get('Title'), loc='left')
    # TODO: check that this is ok. Maybe more customization is needed.
    for axis_name in 'XY':
        getattr(ax, 'set_{}label'.format(axis_name.lower()))(plot_features.get(axis_name + 'Label'), loc='top')

        min_lim = plot_features.get(axis_name + 'Min', min(getattr(h, axis_name.lower() + 'Min')() for h in hist_data))
        max_lim = plot_features.get(axis_name + 'Max', max(getattr(h, axis_name.lower() + 'Max')() for h in hist_data))
        ax.set(**{axis_name.lower() + 'lim': (min_lim, max_lim)})

        # Set log scale
        if plot_features.get('Log' + axis_name):
            ax.set(**{axis_name.lower() + 'scale': 'log'})
            getattr(ax, axis_name.lower() + 'axis').set_major_locator(mpl.ticker.LogLocator(numticks=np.inf))
        
        # Set tick marks frequency
        if (axis_name + 'MajorTickMarks') in plot_features and not plot_features.get('Log' + axis_name):
            base = plot_features.get(axis_name + 'MajorTickMarks')*10**(int(np.log10(max_lim))-1)
            ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base))
        
        if (axis_name + 'MinorTickMarks') in plot_features and not plot_features.get('Log' + axis_name):
            getattr(ax, axis_name.lower() + 'axis').set_minor_locator(mpl.ticker.AutoMinorLocator(
                1+plot_features.get(axis_name + 'MinorTickMarks')))
        
        # Add custom ticks
        if (axis_name + 'CustomMajorTicks') in plot_features:
            ax.set(**{axis_name.lower() + 'xticks': plot_features.get(axis_name + 'CustomMajorTicks')[::2],
                axis_name.lower() + 'ticklabels': plot_features.get('XCustomMajorTicks')[1::2]})
            getattr(ax, 'set_{}ticks'.format(axis_name))([], minor=True)  # Turn off minor ticks
        if (axis_name + 'CustomMinorTicks') in plot_features:
            getattr(ax, 'set_{}ticks'.format(axis_name))(plot_features.get(axis_name + 'CustomMinorTicks'), minor=True)
        if plot_features.get('Plot{}TickLabels'.format(axis_name)) == 0:
            ax.set(**{axis_name.lower() + 'ticklabels': []})

# TODO this might not be useful if all should have separate colorbars
def _create_norm(plot_features, default_zmin, default_zmax):
    zmin = plot_features.get('ZMin', default_zmin)
    zmax = plot_features.get('ZMax', default_zmax)

    # Set log scale
    if plot_features.get('LogZ'):
        norm = mpl.colors.LogNorm(vmin=zmin, vmax=zmax)
    else:
        norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
    return norm


def _create_colorbar(ax, norm, plot_features):
    colorbar_tick_format = mpl.ticker.ScalarFormatter(useMathText=True)
    colorbar_tick_format.set_powerlimits(
        (-plt.rcParams.get('axes.formatter.min_exponent'), plt.rcParams.get('axes.formatter.min_exponent'))
    )
    cbar = ax.get_figure().colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=plt.rcParams['image.cmap']), 
        cax=ax, orientation='vertical', label=plot_features.get('ZLabel'), format=colorbar_tick_format)

    # Set tick marks frequency
    if 'ZMinorTickMarks' in plot_features and not plot_features.get('LogZ'):
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(
            1+plot_features.get('ZMinorTickMarks')))

    if 'ZMajorTickMarks' in plot_features and not plot_features.get('LogZ'):
        base = plot_features.get('ZMajorTickMarks')*10**(int(np.log10(zmax))-1)
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base))

    if plot_features.get('PlotZTickLabels') == 0:
        ax.set_yticklabels([])

    # Add custom ticks
    if 'ZCustomMajorTicks' in plot_features:
        cbar.set_ticks(plot_features.get('ZCustomMajorTicks')[::2])
        cbar.set_ticklabels(plot_features.get('ZCustomMajorTicks')[1::2])
    if plot_features.get('ZCustomMinorTicks') is not None:
        ax.yaxis.set_ticks(plot_features.get('ZCustomMinorTicks'), minor=True)

    return cbar


def _plot_projection(hist_data, fig, hist_features, plot_features):
    plot_ratio = plot_features.get('RatioPlot', True) and len(hist_data) > 1
    ncols = len(hist_data) + 1
    nrows = 2 if plot_ratio else 1
    # TODO Try different values for size of colorbar
    width_ratios = np.array([1] * len(hist_data) + [0.0526])
    # TODO hopefully remove constant
    fig.set_size_inches(np.sum(width_ratios) * plt.rcParams['figure.figsize'][0] * 1.5, plt.rcParams['figure.figsize'][1] * ncols)
    
    gs = fig.add_gridspec(ncols=ncols, nrows=nrows, width_ratios=width_ratios)
    
    zmin, zmax = min(h.zMin() for h in hist_data), max(h.zMax() for h in hist_data)
    normalizer = _create_norm(plot_features, min(h.zMin() for h in hist_data), max(h.zMax() for h in hist_data))
    
    for i, (h, settings) in enumerate(zip(hist_data, hist_features)):
        ax = fig.add_subplot(gs[0, i])
        x_edges, y_edges = np.meshgrid(h.xEdges(), h.yEdges())
        # TODO: should colorbar have min and max for all of them instead? In this case, only one colorbar is used.
        im = ax.pcolormesh(x_edges, y_edges, h.zVals(asgrid=True).T, norm=normalizer)    
        _configure_plot(ax, plot_features, hist_data)
        ax.set_title(settings['Title'])

    _create_colorbar(fig.add_subplot(gs[0, -1]), normalizer, plot_features)

    if plot_features.get('RatioPlot', True):
        ref_h = hist_data[0]
        zmin = zmax = np.nan
        for i, h in enumerate(hist_data[1:]):
            ax = fig.add_subplot(gs[1, i+1])
            x_edges, y_edges = np.meshgrid(h.xEdges(), h.yEdges())
            
            z_ratio = (h.zVals(asgrid=True) / ref_h.zVals(asgrid=True)).T
            z_ratio[~np.isfinite(z_ratio)] = np.nan
            zmin_new, zmax_new = np.nanmin(z_ratio), np.nanmax(z_ratio)
            zmin, zmax = np.nanmin([zmin, zmin_new]), np.nanmax([zmax, zmax_new])            
            im = ax.pcolormesh(x_edges, y_edges, z_ratio)
            # TODO format ratio plot
            #  axis limits cannot be set here since it is unknown
        # TODO: format colorbar
        # TODO move
        colorbar_tick_format = mpl.ticker.ScalarFormatter(useMathText=True)
        colorbar_tick_format.set_powerlimits(
            (-plt.rcParams.get('axes.formatter.min_exponent'), plt.rcParams.get('axes.formatter.min_exponent'))
        )
    
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
