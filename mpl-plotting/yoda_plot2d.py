import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
# The code below is used to suppress ratio_surface warning. Another more clear warning will be included instead. 
import warnings
warnings.filterwarnings('ignore', r'.*Z contains NaN values\. This may result in rendering artifacts\..*')


def _histo2d_to_np(h, xy_type):
    """Convert a yoda 2D histogram into 2D arrays corresponding to x, y, and z.

    Parameters
    ----------
    h : Union[yoda.Scatter3D, yoda.Histo2D, yoda.Profile2D]
        The 2D histogram object from which the edges and z values will be extracted
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
    
    Raises
    ------
    ValueError
        If the value of xy_type is invalid.
    
    Notes
    -----
    TODO code works using meshgrid but is probably slower and not as clean as using reshape
    """
    hs = h.mkScatter()

    nrows = np.argmax(hs.yMins()) + 1

    if xy_type == 'edge':
        y_1d = np.append(hs.yMins()[:nrows], hs.yMax())
        x_1d = np.append(hs.xMins()[::nrows], hs.xMax())
    elif xy_type == 'mid':
        y_1d = hs.yVals()[:nrows]
        x_1d = hs.xVals()[::nrows]
    else:
        raise ValueError('Expected the input parameter of "xy_type" to be "mid" or "edge" but got "{}".'.format(xy_type))
    
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
    **kwargs
        All settings for axis. See docstring in `heatmap` but remove the axis name from the argument.
        E.g., xlabel passed to `heatmap` is equvalent to label in `format_axis`.
    """
    if ax is None:
        ax = plt.gca()
    
    axis_name = axis_name.lower()
    # TODO remove these 2 since they are redundant?
    getattr(ax, 'set_{}label'.format(axis_name))(label)
    if lim is not None:
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
    """Small wrapper function to create a matplotlib norm object that is used as the limits for e.g. a colorbar.
    
    Parameters
    ----------
    zmin, zmax : float
        Minimum and maximum value of the norm
    log : bool
        If True, use a logarithmic norm.
    
    Returns
    -------
    Union[mpl.colors.LogNorm, mpl.colors.Normalize]
        The norm object that has been created.
    """
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
        Additional arguments passed to plt.colorbar. See the documentation of fig.colorbar.

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


def _get_axis_kw(kwargs):
    """Create 3 separate dicts of kwargs that can be passed to format_axis.
    The function also removes the format_axis elements from kwargs. 

    Parameters
    ----------
    kwargs : dict
        Keyword arguments from a plotting function, e.g. `heatmap`.

    Returns
    -------
    all_axis_kw : list[dict[str, Any]]
        A list with dicts of keyword arguments for format_axis.
        The elements correspond to kwargs for the x, y, and z axis respectively.
    """
    all_axis_kws = []
    format_axis_kw_names = {
        'label': '%slabel', 'lim': '%slim', 'log': 'log%s',
        'major_ticks': '%smajor_ticks', 'minor_ticks': '%sminor_ticks', 
        'custom_major_ticks': '%scustom_major_ticks', 'custom_minor_ticks': '%scustom_minor_ticks', 
        'plot_ticklabels': 'plot_%sticklabels'
    }
    for axis_name in 'xyz':
        axis_kws = {}
        for format_axis_name, plot_kw_name in format_axis_kw_names.items():
            kw = plot_kw_name % axis_name
            if kw in kwargs:
                axis_kws[format_axis_name] = kwargs.pop(kw)
        all_axis_kws.append(axis_kws)
    return all_axis_kws


def heatmap(yoda_hist, ax=None, showzero=True, colorbar=True, cmap=None, **kwargs):
    """Plot a 2D heatmap plot of a yoda 2D histogram.

    Parameters
    ----------
    yoda_hist : Union[yoda.Scatter3D, yoda.Histo2D, yoda.Profile2D]
        Yoda 2D histogram that will be plotted.
    ax : matplotlib.axes.Axes
        Axes object in which the 2D heatmap plot will be plotted in. 
        If None, use plt.gca().
    showzero : bool
        If True, plot 0-count bins. If False, do not color them.
    colorbar : bool
        If True, plot a color bar.
    cmap
        Colormap that will be used for the heatmap plot and the color bar.
        Can be of any type that matplotlib accepts.
        If None, the default settings in matplotlib will be used.
    pcm_kw : Optional[dict]
        Additional keyword arguments, as a dict, directly passed to ax.pcolormesh, which is called under the hood.
        If None, nothing passed to the funciton.
    cbar_kw : Optional[dict]
        Additional  keyword arguments, besides cmap and ax, passed to `fig.colorbar`.
        Here, fig is the figure that contains ax.

    Other Parameters
    ----------------
    xlabel, ylabel, zlabel : str, optional
        Axis label.
    xlim, ylim, zlim : tuple, optional
        Lower and upper limits of axis, in that order.
    logx, logy, logz : bool, optional
        If True, set the scale of the axis to log scale. By default False
    xmajor_ticks, ymajor_ticks, zmajor_ticks : int, optional
        Digit of the major ticks, by default None
    xminor_ticks, yminor_ticks, zminor_ticks : int, optional
        Number of minor ticks, by default None
    xcustom_major_ticks, ycustom_major_ticks, zcustom_major_ticks : Sequence, optional
        Location of tick, followed by the corresponding tick label.
    xcustom_minor_ticks, ycustom_minor_ticks, zcustom_minor_ticks : Sequence[float], optional
        A list of the locations of minor ticks at arbitrary positions. 
    plot_xticklabels, plot_yticklabels, plot_zticklabels : bool, optional
        If False, do not plot any tick labels, by default True
    **kwargs : Any
        Additional parameters passed to `ax.set`, e.g., `title`.
    
    Returns
    -------
    im
        The matplotlib artist returned by `ax.pcolormesh`.
    cbar
        The colorbar object created. Only returned if colorbar is True. 
    """
    x_edges, y_edges, z_vals = _histo2d_to_np(yoda_hist, 'edge')

    if not showzero:
        z_vals[z_vals==0] = np.nan
    
    return _heatmap_base(x_edges, y_edges, z_vals, ax=ax, showzero=showzero, colorbar=colorbar, cmap=cmap, **kwargs)


def ratio_heatmap(main_hist, ref_hist, ax=None, showzero=True, colorbar=True, cmap=None, **kwargs):
    """Plot a 2D heatmap plot of the ratio between 2 histograms.
    The ratio is calculated as main_hist / ref_hist.

    Parameters
    ----------
    main_hist : Union[yoda.Scatter3D, yoda.Histo2D, yoda.Profile2D]
        Yoda histogram that will be the dividend.
    ref_hist : Union[yoda.Scatter3D, yoda.Histo2D, yoda.Profile2D]
        Yoda histogram that will be the divisor, i.e. the reference histogram.
    ax : matplotlib.axes.Axes
        Axes object in which the 2D heatmap plot will be plotted in. 
        If None, use plt.gca().
    showzero : bool
        If True, plot 0-count bins. If False, do not color them.
    colorbar : bool
        If True, plot a color bar.
    cmap
        Colormap that will be used for the heatmap plot and the color bar.
        Can be of any type that matplotlib accepts.
        If None, the default settings in matplotlib will be used.
    
    Other Parameters
    ----------------
    The same parameters used in `heatmap`. See its documentation.

    Returns
    -------
    im
        The matplotlib artist returned by `ax.pcolormesh`.
    cbar
        The colorbar object created. Only returned if colorbar is True. 
    """
    
    x_edges, y_edges, main_z = _histo2d_to_np(main_hist, 'edge')
    ref_z = ref_hist.zVals().reshape(main_z.shape)
    # Set delimiter to nan if they are 0 to remove error.
    ref_z[ref_z == 0]  = np.nan
    z_ratio = main_z / ref_z
    
    return _heatmap_base(x_edges, y_edges, z_ratio, ax=ax, showzero=showzero, colorbar=colorbar, cmap=cmap, **kwargs)


def _heatmap_base(x_edges, y_edges, z_vals, ax=None, showzero=True, colorbar=True, cmap=None, 
    pcm_kw=None, cbar_kw=None, **kwargs):
    """The base function used for all heatmap plotting functions.

    Parameters
    ----------
    x_edges, y_edges : array-like
        The x and y coordinates of the edges of the 2D histogram that will be plotted. 
        Have shape (n+1, m+1).
    z_vals : array-like
        Z values that will be plotted in the heatmap plot.
        For a histogram, this corresponds to the counts in each bin.
        Has shape (n, m).
    
    The documentation for the other parameters can be seen in `heatmap`.

    Returns
    -------
    The documentation for the return arguments can be seen in `heatmap`.
    """
    if ax is None:
        ax = plt.gca()
    # Replace None with an empty dict
    pcm_kw, cbar_kw = [({} if kw is None else kw) for kw in (pcm_kw, cbar_kw)]

    norm = _create_norm(*kwargs.get('zlim', (None, None)), kwargs.get('logz'))
    im = ax.pcolormesh(x_edges, y_edges, z_vals, norm=norm, cmap=cmap, **pcm_kw)

    all_axis_kws = _get_axis_kw(kwargs)
    format_axis('x', ax, **all_axis_kws[0])
    format_axis('y', ax, **all_axis_kws[1])

    if colorbar:
        cbar = _add_colorbar(norm, cmap=cmap, ax=ax, **cbar_kw)
        format_axis('y', cbar.ax, **all_axis_kws[2])
        return_args = im, cbar
    else:
        return_args = im
    ax.set(**kwargs)

    return return_args


def surface(yoda_hist, ax=None, showzero=True, cmap=None, **kwargs):
    """Create a 3D plot with a surface representing a histogram.

    Parameters
    ----------
    yoda_hist : Union[yoda.Scatter3D, yoda.Histo2D, yoda.Profile2D]
        Yoda 2D histogram that will be plotted.
    ax : Axes3D object
        Axes object in which the surface plot will be plotted in. Must be a 3D axes object.
        If None, use `plt.gca(projection='3d')`.
    showzero : bool
        If True, plot 0-count bins. If False, do not plot them.
        Surface plots in matplotlib do not officially support masked values.
        Setting this to False might therefore cause unexpected behavior.
    cmap
        Colormap that will be used for the surface plot.
        Can be of any type that matplotlib accepts.
        If None, the default settings in matplotlib will be used.
    elev : Optional[float]
        The elevation angle, i.e. the "camera position" when viewing the plot. Measured in degrees.
        If None, use the default setting in matplotlib.
    azim : Optional[float]
        The azimuthal angle, i.e. the "camera position" when viewing the plot. Measured in degrees.
        If None, use the default setting in matplotlib.
    psurf_kw : Optional[dict]
        Additional keyword arguments, as a dict, directly passed to ax.plot_surface, which is called under the hood.
        If None, nothing is passed to the funciton.
    
    Other Parameters
    ----------------    
    See the `Other Parameters` section in `heatmap`.

    Returns
    -------
    im
        The matplotlib artist returned by `ax.plot_surface`.

    Notes
    -----
    Logarithmic axes do not work for surface plots, due to a bug in matplotlib.
    If this feature is crucial, it is advised to use `heatmap` instead.
    """
    x, y, z = list(_histo2d_to_np(yoda_hist, 'mid'))

    if not showzero:
        z[z==0] = np.nan
    
    return _surface_base(x, y, z, ax, showzero=showzero, cmap=cmap, **kwargs)


def ratio_surface(main_hist, ref_hist, ax=None, showzero=True, cmap=None, **kwargs):
    """Create a 3D plot with a surface representing a ratio between main_hist and ref_hist.

    Parameters
    ----------
    main_hist : Union[yoda.Scatter3D, yoda.Histo2D, yoda.Profile2D]
        Yoda histogram that will be the dividend.
    ref_hist : Union[yoda.Scatter3D, yoda.Histo2D, yoda.Profile2D]
        Yoda histogram that will be the divisor, i.e. the reference histogram.
    ax : Axes3D object
        Axes object in which the surface plot will be plotted in. Must be a 3D axes object.
        If None, use `plt.gca(projection='3d')`.
    showzero : bool
        If True, plot 0-count bins. If False, do not plot them.
        Surface plots in matplotlib do not officially support masked values.
        Setting this to False might therefore cause unexpected behavior.
    cmap
        Colormap that will color the surface.
        Can be of any type that matplotlib accepts.
        If None, the default settings in matplotlib will be used.
    **kwargs
        The same parameters used in `surface`. See its documentation.

    Returns
    -------
    im
        The matplotlib artist returned by `ax.plot_surface`.
    
    Notes
    -----
    Logarithmic axes do not work for surface plots, due to a bug in matplotlib.
    If this feature is crucial, it is advised to use `heatmap` instead.
    """
    x, y, main_z = _histo2d_to_np(main_hist, 'mid')
    ref_z = ref_hist.zVals().reshape(main_z.shape)
    # Set delimiter to nan if they are 0 to remove error.
    ref_z[ref_z == 0]  = np.nan
    z = main_z / ref_z
    return _surface_base(x, y, z, ax, showzero=showzero, cmap=cmap, **kwargs)


def _surface_base(x, y, z, ax=None, showzero=True, cmap=None, elev=None, azim=None, psurf_kw=None, **kwargs):
    """The base function used for all surface plotting functions.

    Parameters
    ----------
    x, y, z : array-like
        The x, y, and z coordinates of the points on the surface that will be plotted. 
        Have shape (n, m).
        
    The documentation for the other parameters can be seen in `surface`.

    Returns
    -------
    The documentation for the return arguments can be seen in `surface`.
    """
    
    if ax is None:
        ax = plt.gca(projection='3d')
    psurf_kw = {} if psurf_kw is None else psurf_kw
    
    if np.any(np.isnan(z)):
        warnings.warn('Plotting ratios of histograms with 0-count bins as surface plots can result in unexpected rendering results or even errors.\n'
            'It is recommended to use it carefully or to plot a heatmap instead.')
    
    all_axis_kws = _get_axis_kw(kwargs)
    im = ax.plot_surface(x, y, z, cmap=cmap, **psurf_kw)
    
    for axis_name, format_axis_kw in zip('xyz', all_axis_kws):
        format_axis(axis_name, ax, **format_axis_kw)

    ax.view_init(elev=elev, azim=azim)
    ax.set(**kwargs)

    return im
