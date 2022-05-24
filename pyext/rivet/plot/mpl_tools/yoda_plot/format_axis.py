import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


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
