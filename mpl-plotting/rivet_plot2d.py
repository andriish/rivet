"""Create a rivet-style plot with 2D histograms."""
import os
import logging

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

import yoda_plot as yp


def _prepare_mpl(yaml_dict, plot_features, style_path):
    """Change rcParams for matplotlib so that plots get a rivet-style look. 

    Parameters
    ----------
    yaml_dict : dict 
        The dictionary created from .yoda and .plot files. Contains all info about rcParams
    plot_features : dict
        The "plot features" section of the yaml_dict.
    style_path : str
        Path to the mplstyle file. See plot_2Dhist for its default value.
    """
    # rc2d is styling to be applied in conjuction with the rivet default style. Might move this to .mplstyle file
    rc2d = {
        'yaxis.labellocation': 'top', 'image.cmap': 'viridis', 'axes.labelpad': 0.7,
        'figure.figsize': (4.5, 4.41), 'figure.subplot.hspace': plt.rcParams['figure.subplot.wspace'],
        'figure.subplot.bottom': 0.085, 'figure.subplot.top': 0.9, 'figure.subplot.left': 0.12462526766, 'figure.subplot.right': 0.89
    }

    plot_style = yaml_dict.get('style', 'default') + '.mplstyle'
    plt.style.use((os.path.join(style_path, plot_style), rc2d, yaml_dict.get('rcParams', {})))
    plt.rcParams['xtick.top'] = plot_features.get('XTwosidedTicks', True)
    plt.rcParams['ytick.right'] = plot_features.get('YTwosidedTicks', True)


def _post_process_fig(fig, filename, fileformats, title, sup_title_kw=None, savefig_kw=None, clf_kw=None):
    """Save, add title to figure and finally clear the figure.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The matplotlib figure object that will be saved 
    filename : str
        The output file name, without the file extension.
    fileformats : Iterable[str]
        See outputfileformats in plot_2Dhist.
    title : str
        The title of the entire figure.
    suptitle_kw, savefig_kw, clf_kw : dict
        Additional kwargs passed to its respective matplotlib functions.
        If None, nothing is passed to the functions.
    """
    if sup_title_kw is None:
        sup_title_kw = {}
    if savefig_kw is None:
        savefig_kw = {}
    if clf_kw is None:
        clf_kw = {}

    fig.suptitle(title, **sup_title_kw)
    for ext in fileformats:
        fig.savefig('{}.{}'.format(filename, ext), **savefig_kw)
    fig.clf(**clf_kw)


def _get_zlim(hist_data, plot_features):
    """
    Calculate zlim based on histograms and plot_features.

    Parameters
    ----------
    hist_data : list
        All histograms
    plot_features : dict

    Returns
    -------
    zmin, zmax : float
        Limits to the z axis.
    
    Notes
    -----
    This code is a modified version of the original ylim-calculating code in make-plots.
    """
    zmax = plot_features.get('ZMax', max(h.zMax() for h in hist_data))
    
    minzmin = plot_features.get('ZMin', min(h.zMin() for h in hist_data))
    if 'ZMin' not in plot_features:
        zmin = plot_features.get('ZMin')
    elif plot_features.get('LogZ', True):
        zmin = minzmin if plot_features.get('FullRange', True) else max(minzmin, 2e-7*zmax)
    elif plot_features.get('ShowZero', True):
        zmin = min(0, minzmin)
    else:
        zmin = minzmin

    return zmin, zmax


def _get_axis_kw(zmin, zmax, plot_features):
    """Create the xaxis_kw, yaxis_kw, zaxis_kw arguments that will be passed to the yoda plotting API.

    Parameters
    ----------
    zmin : float
        Lower limit of the z axis.
    zmax : float
        Upper limit of the z axis.
    plot_features : dict
        Dict containing all settings for the plot. 
        Some of these settings will be passed to xaxis_kw etc.

    Returns
    -------
    axis_kw : list[dict[str, Any]]
        The x, y, and z axis_kw dicts, returned in that order.
    """
    axis_kw = {}
    for a in 'XYZ':
        if a == 'Z':
            lim = (zmin, zmax)
        else:
            lim = (plot_features.get(a+'Min'), plot_features.get(a+'Max'))
        axis_kw.update(
            {
                a.lower()+'label': plot_features.get(a+'Label'), a.lower()+'lim': lim, 'log'+a.lower(): plot_features.get('Log'+a), 
                a.lower()+'major_ticks': plot_features.get(a+'MajorTickMarks'), a.lower()+'minor_ticks': plot_features.get(a+'MinorTickMarks'), 
                a.lower()+'custom_major_ticks': plot_features.get(a+'CustomMajorTicks'), 
                a.lower()+'custom_minor_ticks': plot_features.get(a+'CustomMinorTicks'), 'plot_%sticklabels' % a.lower(): plot_features.get('Plot%sTickLabels' % a)
            }
        )
    return axis_kw


def raise_2dtype_error(wrong_type):
    """Small wrapper to raise the 2DType error for easier refactoring."""
    raise ValueError('Expected "2DType" in the input file to be "heatmap" or "surface" but was "{}".'.format(wrong_type))


def plot_2Dhist(hist_data, hist_features, yaml_dict, filename, style_path='plot_styles/', outputfileformats=('png','pdf')):
    """Plot 2D histograms in hist_data by formatting the figure according to the specified settings. 

    Parameters
    ----------
    hist_data : list[yoda.Scatter3D]
        All histograms that will be plotted.
    hist_features : dict
        Plot settings for each histogram.
        Currently only supports the "Title" setting but more should be added. 
    yaml_dict : dict
        Plot settings for the entire figure.
    filename : str
        The output file name, without the file extension.
    style_path : str
        Path to the directory with the rivet mplstyle file(s).
        TODO this is a temporary argument and will likely be removed once mplstyle files are added to mplconfigdir.
    outputfileformats : Iterable[str]
        All the file formats, e.g., png, pdf, svg, in which the figure(s) will be exported as.
    """
    logging.captureWarnings(True)
    plot_features = yaml_dict['plot features']
    # Set default scale of z axis to log. This is to have a similar behavior to the y axis scale 
    # This if statement must be kept here rather than using plot_features.get('LogZ', True) everywhere.
    # Otherwise it would be complicated to get the LogZ value independently of LogX, LogZ in e.g., _get_axis_kw 
    # TODO change this default value based on whether the bin-value quantiles fit better to an exp or linear model? 
    if 'LogZ' not in plot_features:
        plot_features['LogZ'] = True
    
    _prepare_mpl(yaml_dict, plot_features, style_path)

    zmin, zmax = _get_zlim(hist_data, plot_features)
    axis_kw = _get_axis_kw(zmin, zmax, plot_features)
    
    ratio = plot_features.get('RatioPlot', False) and len(hist_data) > 1
    ratio_zmin = plot_features.get('RatioPlotZMin', 0.5)
    ratio_zmax = plot_features.get('RatioPlotZMax', 1.4999)
    ratio_axis_kw = _get_axis_kw(ratio_zmin, ratio_zmax, plot_features)
    # Ratio plots will never have a log z axis 
    ratio_axis_kw['logz'] = False
    default_camera_angle = (plot_features.get('3DElev'), plot_features.get('3DAzim'))

    # TODO if possible, refactor this entire if-else-statement for less code duplication
    if plot_features.get('2DIndividual', True):
        fig = plt.figure()
        for yoda_hist, hist_settings in zip(hist_data, hist_features):
            if plot_features.get('2DType', 'heatmap') == 'heatmap':
                ax = fig.add_subplot(111)
                yp.heatmap(yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), colorbar=True, cmap=plot_features.get('2DColormap', plt.rcParams['image.cmap']),
                    cbar_kw=dict(fraction=0.075, pad=0.02, aspect=25), **axis_kw)
            elif plot_features.get('2DType') == 'surface':
                ax = fig.add_subplot(111, projection='3d')
                yp.surface(yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), cmap=plot_features.get('2DColormap', plt.rcParams['image.cmap']),
                    elev=default_camera_angle[0], azim=default_camera_angle[1], **axis_kw)
            else:
                raise_2dtype_error(plot_features['2DType'])
    
            _post_process_fig(fig, '{}-{}'.format(filename, hist_settings['Title']), outputfileformats, plot_features.get('Title'))

        if ratio:
            ref_hist = hist_data[0]
            for yoda_hist, hist_settings in zip(hist_data[1:], hist_features[1:]):
                if plot_features.get('2DType', 'heatmap') == 'heatmap':
                    ax = fig.add_subplot(111)
                    yp.ratio_heatmap(ref_hist, yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), colorbar=True, cmap=plot_features.get('RatioPlot2DColormap', 'PRGn'),
                        cbar_kw=dict(fraction=0.075, pad=0.02, aspect=25), **ratio_axis_kw)
                elif plot_features.get('2DType') == 'surface':
                    ax = fig.add_subplot(111, projection='3d')
                    yp.ratio_surface(ref_hist, yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), cmap=plot_features.get('RatioPlot2DColormap', 'PRGn'),
                        elev=plot_features.get('RatioPlot3DElev', default_camera_angle[0]), azim=plot_features.get('RatioPlot3DAzim', default_camera_angle[1]), **ratio_axis_kw)
                else:
                    raise_2dtype_error(plot_features['2DType'])

                _post_process_fig(fig, '{}-{}-ratio'.format(filename, hist_settings['Title']), outputfileformats, plot_features.get('Title'))

    else:
        ncols = len(hist_data)
        nrows = 1 + ratio
        # The last column is intentionally made larger, since it will contain the color bar as well.
        width_ratios = [1] * (len(hist_data) - 1) + [1.1]
        fig = plt.figure(figsize=np.array(plt.rcParams['figure.figsize']) * np.array([sum(width_ratios), nrows]))
        gs = fig.add_gridspec(ncols=ncols, nrows=nrows, width_ratios=width_ratios)

        for i, (yoda_hist, hist_settings) in enumerate(zip(hist_data, hist_features)):
            if plot_features.get('2DType', 'heatmap') == 'heatmap':
                ax = fig.add_subplot(gs[0, i])
                yp.heatmap(yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), colorbar=len(hist_data)-1 == i, cmap=plot_features.get('2DColormap', plt.rcParams['image.cmap']),
                    cbar_kw=dict(fraction=0.075, pad=0.02, aspect=25), **axis_kw)
            elif plot_features.get('2DType') == 'surface':
                ax = fig.add_subplot(gs[0, i], projection='3d')
                yp.surface(yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), cmap=plot_features.get('2DColormap', plt.rcParams['image.cmap']),
                    elev=default_camera_angle[0], azim=default_camera_angle[1], **axis_kw)
            else:
                raise_2dtype_error(plot_features['2DType'])

            ax.set_title(hist_settings['Title'])
        if ratio:
            ref_hist = hist_data[0]
            for i, (yoda_hist, hist_settings) in enumerate(zip(hist_data[1:], hist_features[1:])):
                if plot_features.get('2DType', 'heatmap') == 'heatmap':
                    ax = fig.add_subplot(gs[1, i+1])
                    yp.ratio_heatmap(ref_hist, yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), colorbar=len(hist_data)-1 == i+1, cmap=plot_features.get('RatioPlot2DColormap', 'PRGn'),
                        cbar_kw=dict(fraction=0.075, pad=0.02, aspect=25), **ratio_axis_kw)
                elif plot_features.get('2DType') == 'surface':
                    ax = fig.add_subplot(gs[1, i+1], projection='3d')
                    yp.ratio_surface(ref_hist, yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), cmap=plot_features.get('RatioPlot2DColormap', 'PRGn'),
                    elev=plot_features.get('RatioPlot3DElev', default_camera_angle[0]), azim=plot_features.get('RatioPlot3DAzim', default_camera_angle[1]), **ratio_axis_kw)
                else: 
                    raise_2dtype_error(plot_features['2DType'])
                # TODO Make this label customizable? 
                ax.set_title('{}/{}'.format(hist_settings['Title'], hist_features[0]['Title']))
        _post_process_fig(fig, filename, outputfileformats, plot_features.get('Title'))
    
    plt.close(fig)
