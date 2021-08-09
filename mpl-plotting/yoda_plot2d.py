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
        ax=ax, orientation='vertical', format=colorbar_tick_format, *args, **kwargs)

    return cbar


def proj(yoda_hist, ax=None, zmin=None, zmax=None, logz=False, showzero=True, 
    colorbar=True, cmap=None, pcm_kw=None, cbar_kw=None):
    """Plot a 2D projection plot of a yoda Scatter3D object.

    Parameters
    ----------
    yoda_hist : yoda.Scatter3D
        Yoda scatter that will be plotted.
    
    ax : matplotlib.axes.Axes
        Axes 
    """
    x_edges, y_edges, z_vals = _scatter_to_2d(yoda_hist, 'edge')

    return _proj_base(x_edges, y_edges, z_vals, ax, zmin, zmax, logz, showzero, colorbar, cmap, pcm_kw, cbar_kw)


def ratio_proj(yoda_hist, ref_hist, ax=None, zmin=None, zmax=None, logz=False, showzero=True, 
    colorbar=True, cmap=None, pcm_kw=None, cbar_kw=None):

    x_edges, y_edges, mc_z = _scatter_to_2d(yoda_hist, 'edge')
    ref_z = ref_hist.zVals().reshape(mc_z.shape) # TODO might not work

    z_ratio = mc_z / ref_z
    
    return _proj_base(x_edges, y_edges, z_ratio, ax, zmin, zmax, logz, showzero, colorbar, cmap, pcm_kw, cbar_kw)
    

def _proj_base(x_edges, y_edges, z_vals, ax=None, zmin=None, zmax=None, logz=False, showzero=True, 
    colorbar=True, cmap=None, pcm_kw=None, cbar_kw=None):
    if ax is None:
        ax = plt.gca()
    if pcm_kw is None:
        pcm_kw = {}
    if cbar_kw is None:
        cbar_kw = {}

    norm = _create_norm(zmin, zmax, logz)
    
    if not showzero:
        z_vals[z_vals==0] = np.nan
        
    im = ax.pcolormesh(x_edges, y_edges, z_vals, norm=norm, cmap=cmap, **pcm_kw)

    if colorbar:
        cbar = _add_colorbar(norm, cmap=cmap, ax=ax, **cbar_kw)
        return im, cbar
    
    # TODO keep format_axis inside here. Ask for preference
    
    return im


def surf(yoda_hist, ax=None, zmin=None, zmax=None, showzero=True, cmap=None, psurf_kw=None):
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

    return _surf_base(x, y, z, ax, zmin, zmax, showzero, cmap, psurf_kw)


def ratio_surf(ref_hist, yoda_hist, ax=None, zmin=None, zmax=None, showzero=True, cmap=None, psurf_kw=None):

    x, y, mc_z = _scatter_to_2d(yoda_hist, 'mid')
    ref_z = ref_hist.zVals().reshape(mc_z.shape)
    z = mc_z / ref_z
    return _surf_base(x, y, z, ax, zmin, zmax, showzero, cmap, psurf_kw)

def _surf_base(x, y, z, ax=None, zmin=None, zmax=None, showzero=True, cmap=None, psurf_kw=None):
    if ax is None:
        ax = plt.gca(projection='3d')

    if not showzero:
        z[z==0] = np.nan
    
    im = ax.plot_surface(x, y, z, cmap=cmap)
    
    # If an axis is log/scaled, change axis scale and change axis limits and ticks accordingly.
    # See Notes in docstring.
    for i, (axis_name, lim) in enumerate(zip('XYZ', ((plot_features.get('XMin'), plot_features.get('XMax')), 
                                                     (plot_features.get('YMin'), plot_features.get('YMax')), 
                                                     (zmin, zmax)))):
        if plot_features.get('Log' + axis_name):
            _set_log_axis_3d(axis_name.lower(), ax)
            log_lim = [np.log10(j) if j is not None else j for j in lim]
            format_axis(axis_name, ax, label=preprocess(plot_features.get(axis_name + 'Label')), lim=log_lim, custom_major_ticks=_preprocess_custom_ticks(plot_features.get(axis_name + 'CustomMajorTicks')), custom_minor_ticks=plot_features.get(axis_name+'CustomMinorTicks'), plot_ticklabels=plot_features.get('Plot{}TickLabels'.format(axis_name)))
        else:
            format_axis(axis_name, ax, preprocess(plot_features.get(axis_name + 'Label')), lim=lim, major_ticks=plot_features.get(axis_name + 'MajorTickMarks'), minor_ticks=plot_features.get(axis_name + 'MinorTickMarks'), custom_major_ticks=_preprocess_custom_ticks(plot_features.get(axis_name + 'CustomMajorTicks')), custom_minor_ticks=plot_features.get(axis_name + 'CustomMinorTicks'), plot_ticklabels=plot_features.get('Plot{}TickLabels'.format(axis_name)))

    # TODO rename keywords?
    ax.view_init(elev=plot_features.get('3DElev'), azim=plot_features.get('3DAzim'))

    return im
