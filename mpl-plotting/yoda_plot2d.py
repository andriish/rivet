import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


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
    TODO code works using meshgrid but is probably slower and not as clean as using reshape
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
    If an optional input variable is None, the option will be ignored, unless a different behavior is specified.
    
    Parameters
    ----------
    axis_name : str
        x, y or z. Either upper or lower case.
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
        Location of tick, followed by the corresponding tick label.
    custom_minor_ticks : Sequence[float], optional
        A list of the locations of minor ticks at arbitrary positions. 
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
    else:
        # Set tick marks frequency
        if major_ticks:
            base = major_ticks*10**(int(np.log10(lim[1]))-1)
            getattr(ax, axis_name + 'axis').set_major_locator(mpl.ticker.MultipleLocator(base))
        
        if minor_ticks:
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

    colorbar_tick_format = mpl.ticker.ScalarFormatter(useMathText=True)
    cbar = ax.get_figure().colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax, format=colorbar_tick_format, *args, **kwargs)

    return cbar


def proj(yoda_hist, ax=None, showzero=True, colorbar=True, cmap=None, **kwargs):
    """Plot a 2D projection plot of a yoda Scatter3D object.

    Parameters
    ----------
    yoda_hist : yoda.Scatter3D
        Yoda scatter that will be plotted.
    
    ax : matplotlib.axes.Axes
        Axes 
    """
    x_edges, y_edges, z_vals = _scatter_to_2d(yoda_hist, 'edge')

    return _proj_base(x_edges, y_edges, z_vals, ax=ax, showzero=showzero, colorbar=colorbar, cmap=cmap, **kwargs)


def ratio_proj(yoda_hist, ref_hist, ax=None, showzero=True, colorbar=True, cmap=None, **kwargs):

    x_edges, y_edges, mc_z = _scatter_to_2d(yoda_hist, 'edge')
    ref_z = ref_hist.zVals().reshape(mc_z.shape) # TODO might not work

    z_ratio = mc_z / ref_z
    
    return _proj_base(x_edges, y_edges, z_ratio, ax=ax, showzero=showzero, colorbar=colorbar, cmap=cmap, **kwargs)
    

def _proj_base(x_edges, y_edges, z_vals, ax=None, showzero=True, colorbar=True, cmap=None, 
    pcm_kw=None, cbar_kw=None, xaxis_kw=None, yaxis_kw=None, zaxis_kw=None):
    if ax is None:
        ax = plt.gca()
    pcm_kw, cbar_kw, xaxis_kw, yaxis_kw, zaxis_kw = [{} if kw is None else kw for kw in (pcm_kw, cbar_kw, xaxis_kw, yaxis_kw, zaxis_kw)]

    norm = _create_norm(*zaxis_kw.get('lim', ()), zaxis_kw.get('log'))
    
    if not showzero:
        z_vals[z_vals==0] = np.nan
        
    im = ax.pcolormesh(x_edges, y_edges, z_vals, norm=norm, cmap=cmap, **pcm_kw)

    # TODO ask for preference about keeping format_axis in here
    format_axis('x', ax, **xaxis_kw)
    format_axis('y', ax, **yaxis_kw)
    
    if colorbar:
        cbar = _add_colorbar(norm, cmap=cmap, ax=ax, **cbar_kw)
        format_axis('y', cbar.ax, **zaxis_kw)
        return im, cbar

    return im


def surf(yoda_hist, ax=None, showzero=True, cmap=None, **kwargs):
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
    """
    x, y, z = list(_scatter_to_2d(yoda_hist, 'mid'))

    return _surf_base(x, y, z, ax, showzero=showzero, cmap=cmap, **kwargs)


def ratio_surf(ref_hist, yoda_hist, ax=None, showzero=True, cmap=None, **kwargs):
    x, y, mc_z = _scatter_to_2d(yoda_hist, 'mid')
    ref_z = ref_hist.zVals().reshape(mc_z.shape)
    z = mc_z / ref_z
    return _surf_base(x, y, z, ax, showzero=showzero, cmap=cmap, **kwargs)


def _surf_base(x, y, z, ax=None, showzero=True, cmap=None, psurf_kw=None, xaxis_kw=None, yaxis_kw=None, zaxis_kw=None, elev=None, azim=None):
    if ax is None:
        ax = plt.gca(projection='3d')
    psurf_kw, xaxis_kw, yaxis_kw, zaxis_kw = [{} if kw is None else kw for kw in (psurf_kw, xaxis_kw, yaxis_kw, zaxis_kw)]

    if not showzero:
        z[z==0] = np.nan
    
    im = ax.plot_surface(x, y, z, cmap=cmap, **psurf_kw)
    
    for axis_name, format_axis_kw in zip('xyz', (xaxis_kw, yaxis_kw, zaxis_kw)):
        format_axis(axis_name, ax, **format_axis_kw)

    ax.view_init(elev=elev, azim=azim)

    return im
