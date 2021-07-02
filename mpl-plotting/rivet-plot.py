"""Module creates a rivet-style plot as a pdf."""
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yoda
from ruamel.yaml import YAML


def apply_style(yaml_dicts):
    # TODO: add other styles and perform error checks
    plot_style = yaml_dicts['plot_style'] + '.mplstyle'
    plt.style.use((os.path.join(sys.path[0], plot_style),
                   yaml_dicts['rcParams']))


def parse_yoda_hist(yaml_dicts):
    # TODO: There is probably a more elegant approach
    histograms = []
    for hist_data in yaml_dicts['histograms'].values():
        temp_file = open("temp_file.txt", "w")
        temp_file.write(hist_data)
        temp_file.close()
        histograms.append(list(yoda.readFLAT("temp_file.txt").values())[0])
        os.remove("temp_file.txt")
    return histograms


def rivet_plot(yaml_file):
    """Create plot from yaml file dictionaries. 

    The yaml file contains rcParams for mpl, histogram data, and plot styles.
    """
    # Parse yaml file for rcParams, histogram data, and plot style.
    with open(yaml_file) as file:
        yaml_dicts = YAML(typ='safe').load(file)
    apply_style(yaml_dicts)
    histograms = parse_yoda_hist(yaml_dicts)
    plot_features = yaml_dicts['plot_features']

    if plot_features.get('RatioPlot'):
        fig, (ax, ax_ratio) = plt.subplots(2, 1, sharex=True,
                                           gridspec_kw={'height_ratios': (2, 1)})
    else:
        fig, ax = plt.subplots(1, 1)
    ax.set_xlabel(plot_features.get('XLabel'))
    ax.set_ylabel(plot_features.get('YLabel'), loc='top')
    ax.set_title(plot_features.get('Title'), loc='left')

    XMin = plot_features.get('XMin', min([h.xMin() for h in histograms]))
    XMax = plot_features.get('XMax', max([h.xMax() for h in histograms]))
    ax.set_xlim(XMin, XMax)
    YMin = plot_features.get('YMin', min([h.yMin() for h in histograms]))
    YMax = plot_features.get('YMax', max([h.yMax() for h in histograms]))
    ax.set_ylim(YMin, YMax)

    # Set log scale
    if plot_features.get('LogX'):
        ax.set_xscale('log')
        ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0))
    if plot_features.get('LogY'):
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0))

    # Create useful variables
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    x_points = (histograms[0].xMins() + histograms[0].xMaxs())/2
    x_bins = np.append(histograms[0].xMins(), histograms[0].xMax())
    data_yVals = histograms[0].yVals()

    # Plot data
    ax.hlines(data_yVals, histograms[0].xMins(),
              histograms[0].xMaxs(), 'k')
    ax.plot(x_points, data_yVals, 'ko')
    data_errminus = [err[0] for err in histograms[0].yErrs()]
    data_errplus = [err[1] for err in histograms[0].yErrs()]
    if plot_features.get('ErrorBars'):
        ax.vlines(x_points, (data_yVals - data_errminus),
                  (data_yVals + data_errplus), 'k')

    # Plot mcs
    for i, mc in enumerate(histograms[1:]):
        color = colors[i % len(colors)]
        y_mc = np.insert(mc.yVals(), 0, mc.yVals()[0])
        ax.plot(x_bins, y_mc, color, drawstyle='steps-pre',
                solid_joinstyle='miter')
        if plot_features.get('ErrorBars'):
            mc_errminus = [err[0] for err in mc.yErrs()]
            mc_errplus = [err[1] for err in mc.yErrs()]
            ax.vlines(x_points, (mc.yVals() - mc_errminus),
                      (mc.yVals() + mc_errplus), color)

    # Create ratio plot
    if plot_features.get('RatioPlot'):  # TODO: Check if errorbar/errorbands
        ax_ratio.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        ax_ratio.set_ylabel(plot_features.get('RatioPlotYLabel', 'MC/Data'))
        RatioPlotYMin = plot_features.get('RatioPlotYMin', 0.5)
        RatioPlotYMax = plot_features.get('RatioPlotYMax', 1.4999)
        ax_ratio.set_ylim(RatioPlotYMin, RatioPlotYMax)

        # Plot data
        ax_ratio.hlines(1, XMin, XMax, 'k', zorder=2)
        ax_ratio.plot(x_points, np.ones(len(x_points)), 'ko', zorder=3)
        if plot_features.get('ErrorBands'):
            errbandminus = np.insert((data_yVals - data_errminus)/data_yVals, 0,
                                     ((data_yVals - data_errminus)/data_yVals)[0])
            errbandplus = np.insert((data_yVals + data_errplus)/data_yVals, 0,
                                    ((data_yVals + data_errplus)/data_yVals)[0])
            ax_ratio.fill_between(x_bins, errbandminus, errbandplus,
                                  step='pre', alpha=0.5, zorder=0)
        elif plot_features.get('ErrorBars'):
            ax_ratio.vlines(x_points, (data_yVals - data_errminus)/data_yVals,
                            (data_yVals + data_errplus)/data_yVals, 'k')

        # Plot mcs
        for i, mc in enumerate(histograms[1:]):
            color = colors[i % len(colors)]
            y_ratio = (np.insert(mc.yVals(), 0, mc.yVals()[0])
                       / np.insert(data_yVals, 0, data_yVals[0]))
            ax_ratio.plot(x_bins, y_ratio, color, drawstyle='steps-pre', zorder=1,
                          solid_joinstyle='miter')
            if plot_features.get('ErrorBars'):
                mc_errminus = [err[0] for err in mc.yErrs()]
                mc_errplus = [err[1] for err in mc.yErrs()]
                ax_ratio.vlines(x_points, (mc.yVals() - mc_errminus)/data_yVals,
                                (mc.yVals() + mc_errplus)/data_yVals, color, zorder=1)

    # Legend
    # TODO: Find a better way to implement the custom Data legend graphic

    if plot_features.get('Legend'):
        handles = [AnyObject()]
        labels = ['Data']
        for i, mc in enumerate(histograms[1:]):
            color = colors[i % len(colors)]
            handles.append(mpl.lines.Line2D([], [], color=color))
            labels.append('mc{}'.format(i+1))
        ax.legend(handles, labels, loc='upper left', bbox_to_anchor=[0.5, 0.97],
                  handler_map={AnyObject: AnyObjectHandler()})

    plt.savefig(os.path.join(sys.path[0], 'mpl_d10-x01-y01.pdf'), dpi=500)


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


if __name__ == '__main__':
    rivet_plot('ATLAS_2010_S8918562_d10-x01-y01.yaml')
