"""Module creates a rivet-style plot as a pdf."""
import os
import sys
import io

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yoda
from yamlio import read_yamlfile
from rivet_plot2d import plot_2Dhist


def rivet_plot(yaml_file, plot_name, outputdir='.'):
    """Create plot from the yaml file.

    Parameters
    ----------
    yaml_file : str or dict
        Either the file name of the yaml file or a dictionary object containing the plot data.
    plot_name : str
        Name of the created plot file.
    outputdir : str, optional
        Name of the relative output directory. Default is the current directory.
    """
    # Parse yaml file for histogram data
    if isinstance(yaml_file, str):  # If the file name of the yaml file is passed
        yaml_dicts = read_yamlfile(yaml_file)
        hist_data = _parse_yoda_hist(yaml_dicts)
    else:  # If the dictionary object is passed
        yaml_dicts = yaml_file
        hist_data = [hist_dict['yoda']
                     for hist_dict in yaml_dicts.get('histograms').values()]

    # TODO the first element in hist_data and hist_features should have IsRef==1
    hist_features = [val for val in yaml_dicts['histograms'].values()]
    output_filename = os.path.join(outputdir, plot_name.strip('/'))

    hist_data_scatter = [h.mkScatter() for h in hist_data]
    
    if all(isinstance(h, (yoda.core.Scatter2D, yoda.core.Histo1D, yoda.core.Profile1D)) for h in hist_data):
        # TODO: Refactor code so 1Dhist and 2Dhist calls are similar
        plot_features = yaml_dicts.get('plot features', {})
        fig, axes = _create_plot(yaml_dicts, plot_features, hist_data_scatter)
        _plot_1Dhist(hist_data_scatter, axes, hist_features, plot_features)
        _save_fig(fig, output_filename)

    elif all(isinstance(h, (yoda.Histo2D, yoda.Scatter3D, yoda.Profile2D)) for h in hist_data):
        plot_2Dhist(hist_data_scatter, hist_features, yaml_dicts, output_filename)
    
    else:
        print('Error with Class types:', [type(h) for h in hist_data])
        raise NotImplementedError('Class type cannot be plotted yet')


def _parse_yoda_hist(yaml_dicts):
    """Read yoda string and return yoda object."""
    hist_data = []
    for hist_dict in yaml_dicts['histograms'].values():
        with io.StringIO(hist_dict['yoda']) as file_like:
            hist_data.append(yoda.readYODA(file_like, asdict=False)[0])
    return hist_data


def _create_plot(yaml_dicts, plot_features, hist_data):
    """Create and return Matplotlib fig and axes objects with modified plot features."""
    plot_style = yaml_dicts['style'] + \
        '.mplstyle'  # TODO: Change location of mplstyle file
    if yaml_dicts.get('rcParams'):
        plt.style.use((os.path.join(sys.path[0], plot_style),
                       yaml_dicts['rcParams']))
    else:
        plt.style.use((os.path.join(sys.path[0], plot_style)))

    plt.rcParams['xtick.top'] = plot_features.get('XTwosidedTicks', True)
    plt.rcParams['ytick.right'] = plot_features.get('YTwosidedTicks', True)

    if plot_features.get('RatioPlot'):
        fig, (ax, ax_ratio) = plt.subplots(2, 1, sharex=True,
                                           gridspec_kw={'height_ratios': (2, 1)})
    else:
        fig, ax = plt.subplots(1, 1)

    # Set text labels
    if plot_features.get('RatioPlot'):
        ax_ratio.set_xlabel(plot_features.get('XLabel'))
    else:
        ax.set_xlabel(plot_features.get('XLabel'))
    ax.set_ylabel(plot_features.get('YLabel'), loc='top')
    ax.set_title(plot_features.get('Title'), loc='left')

    # Set plot lims
    XMin = plot_features.get('XMin', min([h.xMin() for h in hist_data]))
    XMax = plot_features.get('XMax', max([h.xMax() for h in hist_data]))
    ax.set_xlim(XMin, XMax)

    # Use maximum y value from all hist datasets
    max_ymax = max([h.yMax() for h in hist_data])
    if plot_features.get('YMax') is not None:
        YMax = plot_features.get('YMax')
    elif plot_features.get('LogY'):
        YMax = 1.7*max_ymax
    else:
        YMax = 1.1*max_ymax

    # Use minimum y value from all hist datasets
    min_ymin = min([h.yMin() for h in hist_data])
    if plot_features.get('YMin') is not None:
        YMin = plot_features.get('YMin')
    elif plot_features.get('LogY'):
        YMin = (min_ymin/1.7 if plot_features.get('FullRange')
                else max(min_ymin/1.7, 2e-7*YMax))
    elif plot_features.get('ShowZero'):
        YMin = 0 if min_ymin > -1e-4 else 1.1*min_ymin
    else:
        YMin = (1.1*min_ymin if min_ymin < -1e-4 else 0 if min_ymin < 1e-4
                else 0.9*min_ymin)

    ax.set_ylim(YMin, YMax)

    # Set log scale and log tick marks frequency
    if plot_features.get('LogX'):
        ax.set_xscale('log')
        ax.xaxis.set_major_locator(mpl.ticker.LogLocator(numticks=np.inf))
    if plot_features.get('LogY'):
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(mpl.ticker.LogLocator(numticks=np.inf))
        ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(
            base=10.0, subs=[i for i in np.arange(0, 1, 0.1)], numticks=np.inf))

    # Set tick marks frequency given the last digit in the tick mark precision for non-log plots
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

    # Add custom ticks for x and y axes
    if plot_features.get('XCustomMajorTicks') is not None:
        ax.set_xticks(plot_features.get('XCustomMajorTicks')[::2])
        ax.set_xticklabels(plot_features.get('XCustomMajorTicks')[1::2])
        ax.set_xticks([], minor=True)  # Turn off minor xticks
        if plot_features.get('RatioPlot'):
            ax_ratio.set_xticks([], minor=True)
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
    if plot_features.get('RatioPlot'):
        return fig, (ax, ax_ratio)
    return fig, (ax,)


def _plot_1Dhist(hist_data, axes, hist_features, plot_features):
    """Plot the 1D historgram data."""
    # Create useful variables
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    x_points = (hist_data[0].xMins() + hist_data[0].xMaxs())/2
    x_bins = np.append(hist_data[0].xMins(), hist_data[0].xMax())
    data_yVals = hist_data[0].yVals()

    # Create local variables for the list of axes
    if plot_features.get('RatioPlot'):
        ax, ax_ratio = axes
    else:
        ax = axes[0]

    # Plot the reference histogram data
    ax.hlines(data_yVals, hist_data[0].xMins(),
              hist_data[0].xMaxs(), 'k')  # Plot reference data as horizontal lines
    ax.plot(x_points, data_yVals, 'ko')  # Plot black dot in the middle of line
    data_errminus = [err[0] for err in hist_data[0].yErrs()]
    data_errplus = [err[1] for err in hist_data[0].yErrs()]
    if hist_features[0].get('ErrorBars', 1):  # Plot ref data error bars by default
        ax.vlines(x_points, (data_yVals - data_errminus),
                  (data_yVals + data_errplus), 'k')

    # Plot the MC histogram data
    for i, mc in enumerate(hist_data[1:]):
        color = colors[i % len(colors)]  # Cycle through colors for MC hists
        y_mc = np.insert(mc.yVals(), 0, mc.yVals()[0])
        ax.plot(x_bins, y_mc, color, drawstyle='steps-pre',
                solid_joinstyle='miter', zorder=5+i)  # Plot MC hist data
        if hist_features[i].get('ErrorBars', 1):  # Plot MC error bars by default
            mc_errminus =  [err[0] for err in mc.yErrs()]
            mc_errplus = [err[1] for err in mc.yErrs()]
            ax.vlines(x_points, (mc.yVals() - mc_errminus),
                      (mc.yVals() + mc_errplus), color, zorder=5+i)

    # Create ratio plot
    if plot_features.get('RatioPlot'):
        # Ratio plot has y range of 0.5 to 1.5
        ax_ratio.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        ax_ratio.set_ylabel(plot_features.get('RatioPlotYLabel', 'MC/Data'))
        RatioPlotYMin = plot_features.get('RatioPlotYMin', 0.5)
        RatioPlotYMax = plot_features.get(
            'RatioPlotYMax', 1.4999)  # Don't plot 1.5
        ax_ratio.set_ylim(RatioPlotYMin, RatioPlotYMax)

        # Plot reference histogram data
        XMin = plot_features.get('XMin', min([h.xMin() for h in hist_data]))
        XMax = plot_features.get('XMax', max([h.xMax() for h in hist_data]))
        ax_ratio.hlines(1, XMin, XMax, 'k', zorder=2)
        ax_ratio.plot(x_points, np.ones(len(x_points)), 'ko', zorder=3)
        if plot_features.get('ErrorBands'):
            errbandminus = np.insert((data_yVals - data_errminus)/data_yVals, 0,
                                     ((data_yVals - data_errminus)/data_yVals)[0])
            errbandplus = np.insert((data_yVals + data_errplus)/data_yVals, 0,
                                    ((data_yVals + data_errplus)/data_yVals)[0])
            ax_ratio.fill_between(x_bins, errbandminus, errbandplus,
                                  step='pre', alpha=0.5, zorder=0)
        elif hist_features[0].get('ErrorBars', 1):
            ax_ratio.vlines(x_points, (data_yVals - data_errminus)/data_yVals,
                            (data_yVals + data_errplus)/data_yVals, 'k')

        # Plot MCs histogram data
        for i, mc in enumerate(hist_data[1:]):
            color = colors[i % len(colors)]
            y_ratio = (np.insert(mc.yVals(), 0, mc.yVals()[0])
                       / np.insert(data_yVals, 0, data_yVals[0]))
            ax_ratio.plot(x_bins, y_ratio, color, drawstyle='steps-pre', zorder=1,
                          solid_joinstyle='miter')
            if hist_features[i].get('ErrorBars', 1):
                ax_ratio.vlines(x_points, (mc.yVals() - mc_errminus)/data_yVals,
                                (mc.yVals() + mc_errplus)/data_yVals, color, zorder=1)

    # Legend
    # TODO: Find a better way to implement the custom Data legend graphic

    if plot_features.get('Legend', 1):
        handles = [AnyObject()]
        labels = [hist_features[0].get('Title', 'Data')]
        for i, mc in enumerate(hist_data[1:]):
            color = colors[i % len(colors)]
            handles.append(mpl.lines.Line2D([], [], color=color))
            labels.append(hist_features[0].get('Title', 'mc{}'.format(i+1)))
        if plot_features.get('LegendAlign') is None or plot_features.get('LegendAlign') == 'l':
            legend_pos = (plot_features.get('LegendXPos', 0.5),
                          plot_features.get('LegendYPos', 0.97))
            ax.legend(handles, labels, loc='upper left', bbox_to_anchor=legend_pos,
                      handler_map={AnyObject: AnyObjectHandler()})
        if plot_features.get('LegendAlign') == 'r':
            legend_pos = (plot_features.get('LegendXPos', 0.97),
                          plot_features.get('LegendYPos', 0.97))
            ax.legend(handles, labels, loc='upper right', bbox_to_anchor=legend_pos,
                      handler_map={AnyObject: AnyObjectHandler()}, markerfirst=False)


def _save_fig(fig, plot_name):
    fig.savefig(plot_name+'.pdf')
    fig.savefig(plot_name+'.png')
    plt.close()


class AnyObject(object):
    """Necessary for custom legend handler."""
    pass


class AnyObjectHandler(object):
    """Creates custom legend handler for Data."""

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        patch = mpl.patches.Circle(
            (x0+width/2, y0+height/2), (2.5)**0.5, facecolor='black')
        handlebox.add_artist(patch)
        patch = mpl.patches.Rectangle(
            (-0.4, 3), width+0.8, 0.8, facecolor='black')
        handlebox.add_artist(patch)
        patch = mpl.patches.Rectangle(
            (width/2-0.4, 0), 0.8, height, facecolor='black')
        handlebox.add_artist(patch)
        return patch
