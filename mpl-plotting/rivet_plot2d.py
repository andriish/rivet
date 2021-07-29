"""Create a rivet-style plot with 2D histograms."""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
# TODO
#   Implement ShowZero
#   Implement contour plot
#   Implement surface plot
#   Option to use different colormap for histograms and ratios (to make it easier to see which one is which) 
# TODO probably remove plot_features as input and replace with individual args once the API has been defined 
#   This will make it easy to move all functions to the yoda plotting API
# TODO docstrings


def format_axis(ax, axis_name, label=None, lim=None, log=False, 
    major_ticks=None, minor_ticks=None, custom_major_ticks=None, custom_minor_ticks=None,
    plot_ticklabels=True):
    """Format an axis (e.g. x axis, **NOT** an Axes object) based on the inputs.
    If an optional input variable is None, the default setting in matplotlib will be used.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The Axes object containing the axis.
    axis_name : str
        x, y or z. Either upper or lower case.
        TODO make this an axis object as input instead?
    label : str, optional
        Axis label.
    label_loc : str
        Location of the label. TODO it is probably not good to add this option, since one has to add many options then. Either allow dicts of kwargs to be passed or change afterwards.
    lim : tuple, optional
        Lower and upper limits of axis, in that order.
    log : bool, optional
        If True, set the scale of the axis to log scale. By default False
    TODO rest of docstring
    major_ticks : [type], optional
        TODO change format of this? Currently assumes the rivet format
    minor_ticks : [type], optional
        [description], by default None
    custom_major_ticks : [type], optional
        [description], by default None
    custom_minor_ticks : [type], optional
        [description], by default None
    plot_ticklabels : bool, optional
        [description], by default True
    """
    axis_name = axis_name.lower()
    # TODO Move loc='top' to somewhere else or make it an input parameter
    getattr(ax, 'set_{}label'.format(axis_name))(label)
    getattr(ax, 'set_{}lim'.format(axis_name))(lim)

    # Set log scale
    if log:
        ax.set(**{axis_name + 'scale': 'log'})
        getattr(ax, axis_name + 'axis').set_major_locator(mpl.ticker.LogLocator(numticks=np.inf))
    
    # Set tick marks frequency
    if major_ticks and not log:
        base = major_ticks*10**(int(np.log10(lim[1]))-1)
        getattr(ax, axis_name + 'axis').set_major_locator(mpl.ticker.MultipleLocator(base))
    
    if minor_ticks and not log:
        getattr(ax, axis_name + 'axis').set_minor_locator(mpl.ticker.AutoMinorLocator(
            1+minor_ticks))
    
    # Add custom ticks
    if custom_major_ticks:
        ax.set(**{axis_name + 'ticks': custom_major_ticks[::2],
            axis_name + 'ticklabels': custom_major_ticks[1::2]})
        getattr(ax, 'set_{}ticks'.format(axis_name))([], minor=True)  # Turn off minor ticks
    if custom_minor_ticks:
        getattr(ax, 'set_{}ticks'.format(axis_name))(custom_minor_ticks, minor=True)
    if plot_ticklabels == 0:
        ax.set(**{axis_name + 'ticklabels': []})


def _create_norm(zmin, zmax, log):
    """Create a mpl norm that is used as the limits for the colorbar and pcolormesh."""

    # Set log scale
    if log:
        return mpl.colors.LogNorm(vmin=zmin, vmax=zmax)
    
    return mpl.colors.Normalize(vmin=zmin, vmax=zmax)


def _add_colorbar(ax, norm, plot_features, *args, **kwargs):
    colorbar_tick_format = mpl.ticker.ScalarFormatter(useMathText=True)
    colorbar_tick_format.set_powerlimits(
        (-plt.rcParams.get('axes.formatter.min_exponent'), plt.rcParams.get('axes.formatter.min_exponent'))
    )
    # TODO change this from rcParams to input arg
    cbar = ax.get_figure().colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=plt.rcParams['image.cmap']), 
        ax=ax, orientation='vertical', label=plot_features.get('ZLabel'), format=colorbar_tick_format, *args, **kwargs)

    return cbar


def _plot_projection(ax, yoda_hist, plot_features, zmin=None, zmax=None):
    """Plot a projection plot using pcolormesh.

    Parameters
    ----------
    ax : matplotlib Axes
        [description]
    yoda_hist : yoda Scatter
        [description]
    plot_features : dict
        TODO remove this and turn into individual kwargs once API has been establited.
    """
    norm = _create_norm(zmin, zmax, plot_features.get('LogZ'))
    
    # TODO xEdges etc might not work with Scatter.
    x_edges, y_edges = np.meshgrid(yoda_hist.xEdges(), yoda_hist.yEdges())
    im = ax.pcolormesh(x_edges, y_edges, yoda_hist.zVals(asgrid=True).T, norm=norm)    

    # TODO probably move out of this function
    format_axis(ax, 'x', plot_features.get('XLabel'), (plot_features.get('XMin'), plot_features.get('XMax')), plot_features.get('LogX'), plot_features.get('XMajorTickMarks'), plot_features.get('XMinorTickMarks'), plot_features.get('XCustomMajorTicks'), plot_features.get('XCustomMinorTicks'), plot_features.get('PlotXTickLabels'))
    format_axis(ax, 'y', plot_features.get('YLabel'), (plot_features.get('YMin'), plot_features.get('YMax')), plot_features.get('LogY'), plot_features.get('YMajorTickMarks'), plot_features.get('YMinorTickMarks'), plot_features.get('YCustomMajorTicks'), plot_features.get('YCustomMinorTicks'), plot_features.get('PlotYTickLabels'))
    
    cbar = _add_colorbar(ax, norm, plot_features, fraction=0.075, pad=0.02, aspect=25)
    format_axis(cbar.ax, 'y', plot_features.get('ZLabel'), (zmin, zmax), plot_features.get('LogZ'), plot_features.get('ZMajorTickMarks'), plot_features.get('ZMinorTickMarks'), plot_features.get('ZCustomMajorTicks'), plot_features.get('ZCustomMinorTicks'), plot_features.get('PlotZTickLabels'))
    
    return im, cbar


def _plot_ratio(ax, ref_hist, yoda_hist, plot_features, zmin=None, zmax=None):
    # TODO code is quite similar to plot_projection. Make into the same function with just data as input? 
    norm = _create_norm(zmin, zmax, log=False)

    x_edges, y_edges = np.meshgrid(yoda_hist.xEdges(), yoda_hist.yEdges())
    z_ratio = (yoda_hist.zVals(asgrid=True) / ref_hist.zVals(asgrid=True)).T
    # TODO depending on ShowZero?
    z_ratio[~np.isfinite(z_ratio)] = np.nan
    im = ax.pcolormesh(x_edges, y_edges, z_ratio, norm=norm)

    # TODO probably move out of this function
    format_axis(ax, 'x', plot_features.get('XLabel'), (plot_features.get('XMin'), plot_features.get('XMax')), plot_features.get('LogX'), plot_features.get('XMajorTickMarks'), plot_features.get('XMinorTickMarks'), plot_features.get('XCustomMajorTicks'), plot_features.get('XCustomMinorTicks'), plot_features.get('PlotXTickLabels'))
    format_axis(ax, 'y', plot_features.get('YLabel'), (plot_features.get('YMin'), plot_features.get('YMax')), plot_features.get('LogY'), plot_features.get('YMajorTickMarks'), plot_features.get('YMinorTickMarks'), plot_features.get('YCustomMajorTicks'), plot_features.get('YCustomMinorTicks'), plot_features.get('PlotYTickLabels'))
    
    # TODO colormap setting here. Must remove _add_colorbar plot_features as input for custom cmap to work.
    cbar = _add_colorbar(ax, norm, plot_features, fraction=0.075, pad=0.02, aspect=25)
    format_axis(cbar.ax, 'y', plot_features.get('RatioPlotZLabel', 'MC/Data'), (zmin, zmax), major_ticks=1, plot_ticklabels=plot_features.get('RatioPlotTickLabels'))
    
    return im, cbar


def _post_process_fig(fig, filename, title, sup_title_kw=None, savefig_kw=None, clf_kw=None):
    """Save close, and add title to figure.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        [description]
    filename : str
        [description]
    title : str
        [description]
    """
    if sup_title_kw is None:
        sup_title_kw = {}
    if savefig_kw is None:
        savefig_kw = {}
    if clf_kw is None:
        clf_kw = {}

    fig.suptitle(title, **sup_title_kw)
    # TODO remove args here
    fig.savefig(filename, **savefig_kw)
    fig.clf(**clf_kw)


def _prepare_mpl(plot_features, style_path):
    # Styling to be applied in conjuction with the rivet default style. Might move this to .mplstyle file
    # TODO axes.labelpad is different for x, y axis.
    #  Maybe directly modify using format_axis, specifically set_xlabel(labelpad=something)?
    rc2d = {
        'yaxis.labellocation': 'top', 'image.cmap': 'jet', 'axes.labelpad': 0.7,
        'figure.figsize': (4.5, 4.41), 'figure.subplot.hspace': plt.rcParams['figure.subplot.wspace'],
        'figure.subplot.bottom': 0.085, 'figure.subplot.top': 0.94, 'figure.subplot.left': 0.12462526766, 'figure.subplot.right': 0.89
    }

    if 'style' in plot_features:
        plot_style = plot_features['style'] + '.mplstyle'
    else:
        plot_style = 'default.mplstyle'
    
    plt.style.use((os.path.join(style_path, plot_style), rc2d, plot_features.get('rcParams', {})))
    
    plt.rcParams['xtick.top'] = plot_features.get('XTwosidedTicks', True)
    plt.rcParams['ytick.right'] = plot_features.get('YTwosidedTicks', True)


def plot_2Dhist(hist_data, hist_features, plot_features, filename=None, individual=True, style_path='.'):
    """Plot 2D histogram in hist_data based 

    Parameters
    ----------
    hist_data : list[yoda.Histo2D | yoda.Profile2D | yoda.Scatter3D] TODO only make it work with Scatter in future?
        All histograms that will be plotted.
    hist_features : dict
        Plot settings for each histogram.
        TODO: currently only supports the "Title" setting but more should be added. 
    plot_features : dict
        Plot settings for the entire axes.
    individual : bool
        TODO make this part of plot_features.
    style_path : str
        Path to the directory with the rivet mplstyle file(s).
        TODO this is a temporary argument and will likely become a constant in rivet instead.
    Returns
    -------
    fig : matplotlib.figure.Figure
        Matplotlib Figure object containing all plots.
    """
    _prepare_mpl(plot_features, style_path)

    # TODO move?
    zmin = plot_features.get('ZMin', min(h.zMin() for h in hist_data))
    zmax = plot_features.get('ZMax', max(h.zMax() for h in hist_data))
    ratio = plot_features.get('RatioPlot', True) and len(hist_data) > 1
    ratio_zmin = plot_features.get('RatioPlotZMin', 0.5)
    ratio_zmax = plot_features.get('RatioPlotZMax', 1.4999)
    # TODO remove/replace with something else after testing
    tmp_savefig_kw = dict(transparent=False)

    # TODO: probably separate functions for both of these cases.
    # TODO individual should be part of plot_features instead.
    if individual:
        fig = plt.figure()

        for yoda_hist, hist_settings in zip(hist_data, hist_features):
            ax = fig.add_subplot(111)
            # TODO add more options here such as surface plot etc.
            #  All plot styles hopefully require the same amount of axes, where parts of an axes will be used as a colorbar.
            _plot_projection(ax, yoda_hist, plot_features, zmin, zmax)
            _post_process_fig(fig, '{}-{}{}'.format(filename[:filename.rindex('.')], hist_settings['Title'], filename[filename.rindex('.'):]), plot_features.get('Title'), savefig_kw=tmp_savefig_kw)
            
        if ratio:
            # TODO more code sharing between projection and ratio? Keep as is for now
            ref_hist = hist_data[0]
            for yoda_hist, hist_settings in zip(hist_data[1:], hist_features[1:]):
                ax = fig.add_subplot(111)
                # TODO some or many styles are not applied to ratio plot.
                # TODO add more options here such as surface plot etc. See comment above.
                _plot_ratio(ax, ref_hist, yoda_hist, plot_features, zmin=ratio_zmin, zmax=ratio_zmax)
                _post_process_fig(fig, '{}-{}-ratio{}'.format(filename[:filename.rindex('.')], hist_settings['Title'], filename[filename.rindex('.'):]), plot_features.get('Title'), savefig_kw=tmp_savefig_kw)

    else:
        ncols = len(hist_data)
        nrows = 1 + ratio
        width_ratios = [1] * (len(hist_data) - 1) + [1.3]

        fig = plt.figure(figsize=np.array(plt.rcParams['figure.figsize']) * np.array([ncols * sum(width_ratios), nrows]))
        # TODO: change size of last one
        gs = fig.add_gridspec(ncols=ncols, nrows=nrows, width_ratios=width_ratios)
        for i, (yoda_hist, hist_settings) in enumerate(zip(hist_data, hist_features)):
            ax = fig.add_subplot(gs[0, i])
            _plot_projection(ax, yoda_hist, hist_settings, plot_features)

        if ratio:
            ref_hist = hist_data[0]
            for i, (yoda_hist, hist_settings) in enumerate(zip(hist_data[1:], hist_features[1:])):
                ax = fig.add_subplot(gs[1, i])
                _plot_ratio(ax, ref_hist, yoda_hist, plot_features)
        _post_process_fig(fig, filename, plot_features.get('Title'))
    
    # TODO this does not return anything useful when individual=True.
    return fig


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
        raise NotImplementedError('Expected 2DType to be "projection", "surface" or "contour" '
            'but got {}'.format(plot_features.get('2DType')))
    
    fig.suptitle(plot_features.get('Title'))
    # _configure_plot(fig)

    return fig
