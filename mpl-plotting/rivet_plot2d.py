"""Create a rivet-style plot with 2D histograms."""
import os

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

import yoda_plot2d as yp
from mathtext_preprocessor import preprocess


def _preprocess_custom_ticks(ticks):
    """Apply the mathtext preprocessor on a custom tick list.

    Parameters
    ----------
    ticks : Iterable
        An iterable of ticks in the format of XCustomMajorTicks, i.e. [<tick value 1>, <tick label 1>, <tick value 2>, <tick value 2>, ...]

    Returns
    -------
    list
        A list of the tick labels, where all labels that are strings have been preprocessed to work in mathtext.
    """
    if ticks is None:
        return None
    return [preprocess(tick) for tick in ticks]


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
    # TODO axes.labelpad is different for x, y axis. Personally, this setting looks better than original rivet style. 
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
    # TODO change/remove constants
    if plot_features.get('ZMax') is not None:
        zmax = plot_features.get('ZMax')
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
    elif plot_features.get('ShowZero', True):
        zmin = 0 if minzmin > -1e-4 else 1.1*minzmin
    else:
        zmin = (1.1*minzmin if minzmin < -1e-4 else 0 if minzmin < 1e-4
                else 0.9*minzmin)

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
    axis_kw = []
    for a in 'XYZ':
        if a == 'Z':
            lim = (zmin, zmax)
        else:
            lim = (plot_features.get(a+'Min'), plot_features.get(a+'Max'))
        axis_kw.append(
            dict(label=preprocess(plot_features.get(a+'Label')), lim=lim, log=plot_features.get('Log'+a), 
                major_ticks=plot_features.get(a+'MajorTickMarks'), minor_ticks=plot_features.get(a+'MinorTickMarks'), 
                custom_major_ticks=_preprocess_custom_ticks(plot_features.get(a+'CustomMajorTicks')), 
                custom_minor_ticks=plot_features.get(a+'CustomMinorTicks'), plot_ticklabels=plot_features.get('Plot%sTickLabels' % a))
        )
    return axis_kw


def plot_2Dhist(hist_data, hist_features, yaml_dict, filename, style_path='.', outputfileformats=('png',)):
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
    outputfileformats : Iterable[str]
        All the file formats, e.g., png, pdf, svg, in which the figure(s) will be exported as.
    Returns
    -------
    fig : matplotlib.figure.Figure
        Matplotlib Figure object containing all plots.
    """
    plot_features = yaml_dict.get('plot features', {})
    _prepare_mpl(yaml_dict, plot_features, style_path)

    zmin, zmax = _get_zlim(hist_data, plot_features)
    axis_kw = _get_axis_kw(zmin, zmax, plot_features)
    
    ratio = plot_features.get('RatioPlot', True) and len(hist_data) > 1
    ratio_zmin = plot_features.get('RatioPlotZMin', 0.5)
    ratio_zmax = plot_features.get('RatioPlotZMax', 1.4999)
    ratio_axis_kw = _get_axis_kw(ratio_zmin, ratio_zmax, plot_features)
    # Ratio plots will not have a log z axis 
    ratio_axis_kw[2]['log'] = False

    # TODO if possible, refactor this entire if-else-statement for less code duplication
    if plot_features.get('2DIndividual', True):
        fig = plt.figure()
        for yoda_hist, hist_settings in zip(hist_data, hist_features):
            if plot_features.get('2DType', 'projection') == 'projection':
                ax = fig.add_subplot(111)
                yp.proj(yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), colorbar=True, cmap=plot_features.get('2DColormap', plt.rcParams['image.cmap']),
                    cbar_kw=dict(fraction=0.075, pad=0.02, aspect=25), xaxis_kw=axis_kw[0], yaxis_kw=axis_kw[1], zaxis_kw=axis_kw[2])
            elif plot_features.get('2DType') == 'surface':
                ax = fig.add_subplot(111, projection='3d')
                yp.surf(yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), cmap=plot_features.get('2DColormap', plt.rcParams['image.cmap']),
                    elev=plot_features.get('3DElev'), azim=plot_features.get('3DAzim'), xaxis_kw=axis_kw[0], yaxis_kw=axis_kw[1], zaxis_kw=axis_kw[2])
            else:
                raise ValueError('Expected "2DType" in the input file to be "projection" or "surface" but was "{}".'.format(plot_features['2DType']))
    
            _post_process_fig(fig, '{}-{}'.format(filename, hist_settings['Title']), outputfileformats, preprocess(plot_features.get('Title')))

        if ratio:
            ref_hist = hist_data[0]
            for yoda_hist, hist_settings in zip(hist_data[1:], hist_features[1:]):
                if plot_features.get('2DType', 'projection') == 'projection':
                    ax = fig.add_subplot(111)
                    yp.ratio_proj(ref_hist, yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), colorbar=True, cmap=plot_features.get('RatioPlot2DColormap'),
                        cbar_kw=dict(fraction=0.075, pad=0.02, aspect=25), xaxis_kw=ratio_axis_kw[0], yaxis_kw=ratio_axis_kw[1], zaxis_kw=ratio_axis_kw[2])
                elif plot_features.get('2DType') == 'surface':
                    ax = fig.add_subplot(111, projection='3d')
                    yp.ratio_surf(ref_hist, yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), cmap=plot_features.get('RatioPlot2DColormap', plt.rcParams['image.cmap']),
                        elev=plot_features.get('3DRatioPlotElev'), azim=plot_features.get('RatioPlot3DAzim'), xaxis_kw=ratio_axis_kw[0], yaxis_kw=ratio_axis_kw[1], zaxis_kw=ratio_axis_kw[2])
                _post_process_fig(fig, '{}-{}-ratio'.format(filename, hist_settings['Title']), outputfileformats, preprocess(plot_features.get('Title')))

    else:
        ncols = len(hist_data)
        nrows = 1 + ratio
        # The last column is intentionally made larger, since it will contain the color bar as well.
        width_ratios = [1] * (len(hist_data) - 1) + [1.1]
        fig = plt.figure(figsize=np.array(plt.rcParams['figure.figsize']) * np.array([sum(width_ratios), nrows]))
        gs = fig.add_gridspec(ncols=ncols, nrows=nrows, width_ratios=width_ratios)

        for i, (yoda_hist, hist_settings) in enumerate(zip(hist_data, hist_features)):
            if plot_features.get('2DType', 'projection') == 'projection':
                ax = fig.add_subplot(gs[0, i])
                yp.proj(yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), colorbar=len(hist_data)-1 == i, cmap=plot_features.get('2DColormap', plt.rcParams['image.cmap']),
                    cbar_kw=dict(fraction=0.075, pad=0.02, aspect=25), xaxis_kw=axis_kw[0], yaxis_kw=axis_kw[1], zaxis_kw=axis_kw[2])
            elif plot_features.get('2DType') == 'surface':
                ax = fig.add_subplot(gs[0, i], projection='3d')
                yp.surf(yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), cmap=plot_features.get('2DColormap', plt.rcParams['image.cmap']),
                    elev=plot_features.get('3DElev'), azim=plot_features.get('3DAzim'), xaxis_kw=axis_kw[0], yaxis_kw=axis_kw[1], zaxis_kw=axis_kw[2])
            else:
                raise ValueError('Expected "2DType" in the input file to be "projection" or "surface" but was "{}".'.format(plot_features['2DType']))

            ax.set_title(preprocess(hist_settings['Title']))
        if ratio:
            ref_hist = hist_data[0]
            for i, (yoda_hist, hist_settings) in enumerate(zip(hist_data[1:], hist_features[1:])):
                if plot_features.get('2DType', 'projection') == 'projection':
                    ax = fig.add_subplot(gs[1, i+1])
                    yp.ratio_proj(ref_hist, yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), colorbar=len(hist_data)-1 == i+1, cmap=plot_features.get('RatioPlot2DColormap', plt.rcParams['image.cmap']),
                        cbar_kw=dict(fraction=0.075, pad=0.02, aspect=25), xaxis_kw=ratio_axis_kw[0], yaxis_kw=ratio_axis_kw[1], zaxis_kw=ratio_axis_kw[2])
                elif plot_features.get('2DType') == 'surface':
                    ax = fig.add_subplot(gs[1, i+1], projection='3d')
                    yp.ratio_surf(ref_hist, yoda_hist, ax=ax, showzero=plot_features.get('ShowZero', True), cmap=plot_features.get('RatioPlot2DColormap', plt.rcParams['image.cmap']),
                    elev=plot_features.get('RatioPlot3DElev'), azim=plot_features.get('RatioPlot3DAzim'), xaxis_kw=ratio_axis_kw[0], yaxis_kw=ratio_axis_kw[1], zaxis_kw=ratio_axis_kw[2])
                # TODO Make this label customizable? 
                ax.set_title(preprocess('{}/{}'.format(hist_settings['Title'], hist_features[0]['Title'])))
        _post_process_fig(fig, filename, outputfileformats, preprocess(plot_features.get('Title')))
    
    plt.close(fig)
