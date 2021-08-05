"""Create a rivet-style plot with 2D histograms."""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
# TODO
#   Use mathtext preprocessor on all labels. Some options for where this should be done:
#   - Modify the dict directly once it has been read
#   - Change labels directly inside the plotting functions
#   - Change labels right before passing them into the plotting functions
#   - Change labels in rivet_makeyaml
#   Refactor
#   Probably remove plot_features as input arg at many places and replace with individual args once the API has been defined
#       This will make it easy to move all functions to the yoda plotting API
#   docstrings

def _scatter_to_2d(hs, xy_type):
    """Convert a yoda.Scatter3D into 2D arrays corresponding to x, y, and z.

    Parameters
    ----------
    hs : yoda.Scatter3D
        The scatter object from which the edges z values will be extracted
    xy_type : str
        Either "edge" or "mid". If "edge", return the edges of the of the original bins in the histogram. 
        If "mid", return the midpoints of the x and y bins.
    
    Returns
    -------
    x, y : 2D np.ndarray
        The edges or midpoints of the original histogram, as a 2D array, similar to the output of np.meshgrid. 
        If xy_type is "edge", the shapes of x and y are (n+1, m+1). If xy_type is "mid", the shapes are (n, m).
        Here, n is the number of x bins and m is the number of y bins in the histogram.
    z : 2D np.ndarray
        The z values in each bin of the original histogram, as a 2D array.
        Shape: (n, m), where n is the number of x bins and m is the number of y bins in the histogram.
    
    Notes
    -----
    TODO code works using meshgrid but probably slower and not as clean as using reshape
    """
    nrows = np.argmax(hs.yMins()) + 1
    if xy_type == 'edge':
        y_1d = np.append(hs.yMins()[:nrows], hs.yMax())
        x_1d = np.append(hs.xMins()[::nrows], hs.xMax())
    elif xy_type == 'mid':
        y_1d = hs.yVals()[:nrows]
        x_1d = hs.xVals()[::nrows]
    y, x = np.meshgrid(y_1d, x_1d)
    z = hs.zVals().reshape((-1, nrows))
    return x, y, z


def format_axis(axis_name, ax=None, label=None, lim=None, log=False, 
    major_ticks=None, minor_ticks=None, custom_major_ticks=None, custom_minor_ticks=None,
    plot_ticklabels=True):
    """Format an axis (e.g. x axis, **NOT** an Axes object) based on the inputs.
    If an optional input variable is None, the default setting in matplotlib will be used.
    
    Parameters
    ----------
    axis_name : str
        x, y or z. Either upper or lower case.
        TODO make this an axis object as input instead?
    ax : matplotlib.axes.Axes
        The Axes object containing the axis. If None, uses plt.gca()
    label : str, optional
        Axis label.
    lim : tuple, optional
        Lower and upper limits of axis, in that order.
    log : bool, optional
        If True, set the scale of the axis to log scale. By default False
    major_ticks : int, optional
        Digit of the major ticks, by default None
    minor_ticks : int, optional
        Number of minor ticks, by default None
    custom_major_ticks : Sequence, optional
        [description], by default None
    custom_minor_ticks : Sequence[float], optional
        [description], by default None
    plot_ticklabels : bool, optional
        If False, do not plot any tick labels, by default True
    """
    if ax is None:
        ax = plt.gca()
    
    axis_name = axis_name.lower()
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
    """Small wrapper function to create mpl norm that is used as the limits for the colorbar and pcolormesh."""
    if log:
        return mpl.colors.LogNorm(vmin=zmin, vmax=zmax)
    return mpl.colors.Normalize(vmin=zmin, vmax=zmax)


def _add_colorbar(norm, cmap, ax=None, *args, **kwargs):
    """Small wrapper function to create a colorbar

    Parameters
    ----------
    norm : norm
        Contains information about min, max and scale (log, lin etc) of colorbar.
    cmap : str
        Colormap for the colorbar
    ax : mpl.axes.Axes
        Axes object that the matplotlib will take a part of to create a colorbar.
        If None, uses plt.gca().
    args, kwargs
        Additional arguments passed to plt.colorbar. See its documentation.

    Returns
    -------
    Colorbar
    """
    if ax is None:
        ax = plt.gca()

    # BUG exponent (i.e. 10^x) placement is weird
    colorbar_tick_format = mpl.ticker.ScalarFormatter(useMathText=True)
    colorbar_tick_format.set_powerlimits(
        (-plt.rcParams.get('axes.formatter.min_exponent'), plt.rcParams.get('axes.formatter.min_exponent'))
    )
    cbar = ax.get_figure().colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax, orientation='vertical', format=colorbar_tick_format, *args, **kwargs)

    return cbar


def _plot_projection(yoda_hist, plot_features, ax=None, zmin=None, zmax=None, colorbar=True):
    """Plot a projection plot using pcolormesh.

    Parameters
    ----------
    yoda_hist : yoda.Scatter3D
        Yoda scatter that will be plotted.
    plot_features : dict
        [description]
    ax : matplotlib.axes.Axes
        [description]
    """
    if ax is None:
        ax = plt.gca()

    norm = _create_norm(zmin, zmax, plot_features.get('LogZ', False))
    
    x_edges, y_edges, z_vals = _scatter_to_2d(yoda_hist, 'edge')
    
    if not plot_features.get('ShowZero', True):
        z_vals[z_vals==0] = np.nan
        
    im = ax.pcolormesh(x_edges, y_edges, z_vals, norm=norm, cmap=plot_features.get('2DColormap', 'jet'))

    # TODO probably move out of this function
    format_axis('x', ax, plot_features.get('XLabel'), (plot_features.get('XMin'), plot_features.get('XMax')), plot_features.get('LogX'), plot_features.get('XMajorTickMarks'), plot_features.get('XMinorTickMarks'), plot_features.get('XCustomMajorTicks'), plot_features.get('XCustomMinorTicks'), plot_features.get('PlotXTickLabels'))
    format_axis('y', ax, plot_features.get('YLabel'), (plot_features.get('YMin'), plot_features.get('YMax')), plot_features.get('LogY'), plot_features.get('YMajorTickMarks'), plot_features.get('YMinorTickMarks'), plot_features.get('YCustomMajorTicks'), plot_features.get('YCustomMinorTicks'), plot_features.get('PlotYTickLabels'))
    
    if colorbar:
        cbar = _add_colorbar(norm, cmap=plot_features.get('2DColormap', 'jet'), ax=ax, fraction=0.075, pad=0.02, aspect=25)
        format_axis('y', cbar.ax, plot_features.get('ZLabel'), (zmin, zmax), plot_features.get('LogZ', False), plot_features.get('ZMajorTickMarks'), plot_features.get('ZMinorTickMarks'), plot_features.get('ZCustomMajorTicks'), plot_features.get('ZCustomMinorTicks'), plot_features.get('PlotZTickLabels'))
        return im, cbar
    
    return im


def _plot_ratio_projection(ref_hist, yoda_hist, plot_features, ax=None, zmin=None, zmax=None, colorbar=True):
    # TODO code is quite similar to plot_projection. Make into the same function with just data as input? Keep as is for now
    if ax is None:
        ax = plt.gca()

    norm = _create_norm(zmin, zmax, log=False)

    x_edges, y_edges, mc_z = _scatter_to_2d(yoda_hist, 'edge')
    _, _, ref_z = _scatter_to_2d(ref_hist, 'edge')
    z_ratio = mc_z / ref_z

    if not plot_features.get('ShowZero', True):
        z_ratio[z_ratio==0] = np.nan

    im = ax.pcolormesh(x_edges, y_edges, z_ratio, norm=norm, cmap=plot_features.get('2DRatioColormap', 'jet'))

    # TODO probably move out of this function
    format_axis('x', ax, plot_features.get('XLabel'), (plot_features.get('XMin'), plot_features.get('XMax')), plot_features.get('LogX'), plot_features.get('XMajorTickMarks'), plot_features.get('XMinorTickMarks'), plot_features.get('XCustomMajorTicks'), plot_features.get('XCustomMinorTicks'), plot_features.get('PlotXTickLabels'))
    format_axis('y', ax, plot_features.get('YLabel'), (plot_features.get('YMin'), plot_features.get('YMax')), plot_features.get('LogY'), plot_features.get('YMajorTickMarks'), plot_features.get('YMinorTickMarks'), plot_features.get('YCustomMajorTicks'), plot_features.get('YCustomMinorTicks'), plot_features.get('PlotYTickLabels'))
    
    if colorbar:
        # TODO change default cmap from jet to something else so that ratio plots look different?
        cbar = _add_colorbar(norm, cmap=plot_features.get('2DRatioColormap', 'jet'), ax=ax, fraction=0.075, pad=0.02, aspect=25)
        format_axis('y', cbar.ax, plot_features.get('RatioPlotZLabel', 'MC/Data'), (zmin, zmax), major_ticks=1, plot_ticklabels=plot_features.get('RatioPlotTickLabels'))
        return im, cbar
    
    return im


def _plot_surface(yoda_hist, plot_features, ax=None, zmin=None, zmax=None, *args, **kwargs):
    """Plot 3D plot with a surface representing a histogram.

    Parameters
    ----------
    yoda_hist : yoda.Scatter3D
        Yoda scatter that will be plotted.
    plot_features : dict
        Settings that will be applied to the entire figure or each axes.
    ax : mpl_toolkits.mplot3d.Axes3D
        The axes in which the surface plot will be plotted in. If None, use the latest used 3D axes.
        Later, I might move colorbar function out of _plot_*_projection and this might not be needed.
    args, kwargs
        Will be ignored. Are only there for compatibility reasons.
    """
    # TODO code is quite similar to plot_surface.
    if ax is None:
        ax = plt.gca(projection='3d')

    x_mids, y_mids, z_vals = _scatter_to_2d(yoda_hist, 'mid')

    if not plot_features.get('ShowZero', True):
        z_vals[z_vals==0] = np.nan
        
    im = ax.plot_surface(x_mids, y_mids, z_vals, cmap=plot_features.get('2DColormap', 'jet'))

    # TODO probably move out of this function
    format_axis('x', ax, plot_features.get('XLabel'), (plot_features.get('XMin'), plot_features.get('XMax')), plot_features.get('LogX'), plot_features.get('XMajorTickMarks'), plot_features.get('XMinorTickMarks'), plot_features.get('XCustomMajorTicks'), plot_features.get('XCustomMinorTicks'), plot_features.get('PlotXTickLabels'))
    format_axis('y', ax, plot_features.get('YLabel'), (plot_features.get('YMin'), plot_features.get('YMax')), plot_features.get('LogY'), plot_features.get('YMajorTickMarks'), plot_features.get('YMinorTickMarks'), plot_features.get('YCustomMajorTicks'), plot_features.get('YCustomMinorTicks'), plot_features.get('PlotYTickLabels'))
    format_axis('z', ax, plot_features.get('ZLabel'), (zmin, zmax), plot_features.get('LogZ', False), plot_features.get('ZMajorTickMarks'), plot_features.get('ZMinorTickMarks'), plot_features.get('ZCustomMajorTicks'), plot_features.get('ZCustomMinorTicks'), plot_features.get('PlotZTickLabels'))
    
    # TODO rename keywords?
    ax.view_init(elev=plot_features.get('3DElev'), azim=plot_features.get('3DAzim'))

    return im


def _plot_ratio_surface(ref_hist, yoda_hist, plot_features, ax=None, zmin=None, zmax=None, *args, **kwargs):
    # args, kwargs will be ignored and are there for compatibility reasons. See _plot_surface
    # TODO code is quite similar to plot_surface.
    if ax is None:
        ax = plt.gca(projection='3d')

    x_mids, y_mids, mc_z = _scatter_to_2d(yoda_hist, 'mid')
    _, _, ref_z = _scatter_to_2d(ref_hist, 'mid')
    z_ratio = mc_z / ref_z
    
    if not plot_features.get('ShowZero', True):
        z_ratio[z_ratio==0] = np.nan
    
    im = ax.plot_surface(x_mids, y_mids, z_ratio, cmap=plot_features.get('2DRatioColormap', 'jet'))

    # TODO probably move out of this function
    format_axis('x', ax, plot_features.get('XLabel'), (plot_features.get('XMin'), plot_features.get('XMax')), plot_features.get('LogX'), plot_features.get('XMajorTickMarks'), plot_features.get('XMinorTickMarks'), plot_features.get('XCustomMajorTicks'), plot_features.get('XCustomMinorTicks'), plot_features.get('PlotXTickLabels'))
    format_axis('y', ax, plot_features.get('YLabel'), (plot_features.get('YMin'), plot_features.get('YMax')), plot_features.get('LogY'), plot_features.get('YMajorTickMarks'), plot_features.get('YMinorTickMarks'), plot_features.get('YCustomMajorTicks'), plot_features.get('YCustomMinorTicks'), plot_features.get('PlotYTickLabels'))
    format_axis('z', ax, plot_features.get('RatioPlotZLabel', 'MC/Data'), (zmin, zmax), major_ticks=1, plot_ticklabels=plot_features.get('RatioPlotTickLabels'))
    
    # TODO remove these parameters and use the same as for main plots?
    ax.view_init(elev=plot_features.get('RatioPlot3DElev'), azim=plot_features.get('RatioPlot3DAzim'))

    return im


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
    fig.savefig(filename, **savefig_kw)
    fig.clf(**clf_kw)


def _prepare_mpl(yaml_dict, plot_features, style_path):
    # Styling to be applied in conjuction with the rivet default style. Might move this to .mplstyle file
    # axes.labelpad is different for x, y axis. Personally, this setting looks better than original rivet style. 
    #   Maybe directly modify using format_axis, specifically set_xlabel(labelpad=something)?
    rc2d = {
        'yaxis.labellocation': 'top', 'image.cmap': 'jet', 'axes.labelpad': 0.7,
        'figure.figsize': (4.5, 4.41), 'figure.subplot.hspace': plt.rcParams['figure.subplot.wspace'],
        'figure.subplot.bottom': 0.085, 'figure.subplot.top': 0.94, 'figure.subplot.left': 0.12462526766, 'figure.subplot.right': 0.89
    }

    plot_style = yaml_dict.get('style', 'default') + '.mplstyle'
    plt.style.use((os.path.join(style_path, plot_style), rc2d, yaml_dict.get('rcParams', {})))
    plt.rcParams['xtick.top'] = plot_features.get('XTwosidedTicks', True)
    plt.rcParams['ytick.right'] = plot_features.get('YTwosidedTicks', True)


def _get_zlim(hist_data, plot_features):
    """
    Calculate zlim based on histograms and plot_features.

    Parameters
    ----------
    hist_data : list
        All histograms
    plot_features : dict
        TODO replace with individual input args later.

    Returns
    -------
    zmin, zmax : float
        Limits to the z axis.
    """
    if plot_features.get('ZMax') is not None:
        zmax = plot_features.get('ZMax')
    # TODO why do the constants 1.7, 1.1, 2e-7, 1e-4, 0.9 exist?
    elif plot_features.get('LogZ'):
        zmax = 1.7*max(h.zMax() for h in hist_data)
    else:
        zmax = 1.1*max(h.zMax() for h in hist_data)

    minzmin = min(h.zMin() for h in hist_data)
    if plot_features.get('ZMin') is not None:
        zmin = plot_features.get('ZMin')
    elif plot_features.get('LogZ'):
        zmin = (minzmin/1.7 if plot_features.get('FullRange')
                else max(minzmin/1.7, 2e-7*zmax))
    elif plot_features.get('ShowZero'):
        zmin = 0 if minzmin > -1e-4 else 1.1*minzmin
    else:
        zmin = (1.1*minzmin if minzmin < -1e-4 else 0 if minzmin < 1e-4
                else 0.9*minzmin)
    
    return zmin, zmax


def plot_2Dhist(hist_data, hist_features, yaml_dict, filename, style_path='.', outputfileformat='png'):
    """Plot 2D histogram in hist_data based 

    Parameters
    ----------
    hist_data : list[yoda.Scatter3D]
        All histograms that will be plotted.
    hist_features : dict
        Plot settings for each histogram.
        Currently only supports the "Title" setting but more should be added. 
    yaml_dict : dict
        Plot settings for the entire figure.
    style_path : str
        Path to the directory with the rivet mplstyle file(s).
        TODO this is a temporary argument and will likely become a constant in rivet instead.
    Returns
    -------
    fig : matplotlib.figure.Figure
        Matplotlib Figure object containing all plots.
    """
    plot_features = yaml_dict.get('plot features', {})
    _prepare_mpl(yaml_dict, plot_features, style_path)

    zmin, zmax = _get_zlim(hist_data, plot_features)
    
    ratio = plot_features.get('RatioPlot', True) and len(hist_data) > 1
    ratio_zmin = plot_features.get('RatioPlotZMin', 0.5)
    ratio_zmax = plot_features.get('RatioPlotZMax', 1.4999)

    if plot_features.get('2DType', 'projection') == 'projection':
        plot_function = _plot_projection
        ratio_function = _plot_ratio_projection
        subplot_kw = dict()
    elif plot_features.get('2DType') == 'surface':
        plot_function = _plot_surface
        ratio_function = _plot_ratio_surface
        subplot_kw = dict(projection='3d')
    else:
        raise ValueError('Expected "2DType" in the input file to be "projection" or "surface" but was {}'.format(plot_features['2DType']))
    
    # TODO: probably separate functions for both of these cases.
    if plot_features.get('2DIndividual', True):
        fig = plt.figure()

        for yoda_hist, hist_settings in zip(hist_data, hist_features):
            ax = fig.add_subplot(111, **subplot_kw)
            plot_function(yoda_hist, plot_features, ax=ax, zmin=zmin, zmax=zmax)
            _post_process_fig(fig, '{}-{}.{}'.format(filename, hist_settings['Title'], outputfileformat), plot_features.get('Title'))
            
        if ratio:
            ref_hist = hist_data[0]
            for yoda_hist, hist_settings in zip(hist_data[1:], hist_features[1:]):
                ax = fig.add_subplot(111, **subplot_kw)
                ratio_function(ref_hist, yoda_hist, plot_features, ax=ax, zmin=ratio_zmin, zmax=ratio_zmax)
                _post_process_fig(fig, '{}-{}-ratio.{}'.format(filename, hist_settings['Title'], outputfileformat), plot_features.get('Title'))

    else:
        ncols = len(hist_data)
        nrows = 1 + ratio
        # The last column is intentionally made 10% larger, since it will contain the color bar as well.
        width_ratios = [1] * (len(hist_data) - 1) + [1.1]
        fig = plt.figure(figsize=np.array(plt.rcParams['figure.figsize']) * np.array([sum(width_ratios), nrows]))
        gs = fig.add_gridspec(ncols=ncols, nrows=nrows, width_ratios=width_ratios)

        for i, (yoda_hist, hist_settings) in enumerate(zip(hist_data, hist_features)):
            ax = fig.add_subplot(gs[0, i], **subplot_kw)
            plot_function(yoda_hist, plot_features, ax=ax, zmin=zmin, zmax=zmax, colorbar=len(hist_data)-1 == i)

        if ratio:
            ref_hist = hist_data[0]
            for i, (yoda_hist, hist_settings) in enumerate(zip(hist_data[1:], hist_features[1:])):
                ax = fig.add_subplot(gs[1, i+1], **subplot_kw)
                ratio_function(ref_hist, yoda_hist, plot_features, ax=ax, zmin=ratio_zmin, zmax=ratio_zmax, colorbar=len(hist_data)-2 == i)
        _post_process_fig(fig, '{}.{}'.format(filename, outputfileformat), plot_features.get('Title'))
    
    # TODO this does not return anything useful when 2DIndividual==True. Keep as is anyway?
    return fig
