import os
import io

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yoda
import yoda_plot1d

def plot_1Dhist(hist_data, hist_features, yaml_dicts, filename):
    """Plot the Scatter1D and Scatter2D objects using Rivet styles.

    Parameters
    ----------
    hist_data : list[yoda.Scatter2D]
        All histograms that will be plotted.
    hist_features : dict
        Plot settings for each histogram.
        Currently only supports the "Title" setting but more should be added.
    yaml_dicts : dict
        Plot settings for the entire figure.
    filename : str
        Name of the saved plot file.
    """
    plot_style = os.path.join('plot_styles', yaml_dicts['style'] + '.mplstyle')
    if not os.path.isfile(plot_style):
        raise NotImplementedError('Plot style file not found.')
    if yaml_dicts.get('rcParams'):  # Apply rcParams to mpl
        plt.style.use((plot_style, yaml_dicts.get('rcParams')))
    else:
        plt.style.use(plot_style)

    plot_features = yaml_dicts.get('plot features', {})
    yoda_type = 'hist' if isinstance(hist_data[0], yoda.Scatter2D) else 'scatter'

    # Create fig and axes
    if plot_features.get('RatioPlot', 1) and yoda_type == 'hist':
        fig, (ax, ax_ratio) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': (2, 1)})
    else:
        fig, ax = plt.subplots(1, 1)

    plot_errorbars = [h.get('ErrorBars', 1) for h in hist_features]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # Plot histogram and ratio using Yoda function
    if yoda_type == 'scatter':
        yoda_plot1d.plot_scatter1D(hist_data, ax=ax)
    else:
        yoda_plot1d.plot_hist(hist_data, True, ax, plot_errorbars, colors)
        if plot_features.get('RatioPlot', 1):
            yoda_plot1d.plot_ratio(hist_data, ax_ratio, plot_errorbars, plot_features.get('ErrorBands'), colors)

    plt.rcParams['xtick.top'] = plot_features.get('XTwosidedTicks', True)
    plt.rcParams['ytick.right'] = plot_features.get('YTwosidedTicks', True)

    # Set text labels
    if plot_features.get('RatioPlot', 1) and yoda_type == 'hist':
        ax_ratio.set_xlabel(plot_features.get('XLabel'))
    else:
        ax.set_xlabel(plot_features.get('XLabel'))
    ax.set_ylabel(plot_features.get('YLabel'), loc='top')
    ax.set_title(plot_features.get('Title'), loc='left')

    # Set plot lims
    if yoda_type == 'hist':
        XMin = plot_features.get('XMin', min([h.xMin() for h in hist_data]))
        XMax = plot_features.get('XMax', max([h.xMax() for h in hist_data]))
        ax.set_xlim(XMin, XMax)

    if plot_features.get('RatioPlot', 1) and yoda_type == 'hist':
        # Ratio plot has y range of 0.5 to 1.5
        ax_ratio.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        ax_ratio.set_ylabel(plot_features.get('RatioPlotYLabel', 'MC/Data'))
        RatioPlotYMin = plot_features.get('RatioPlotYMin', 0.5)
        RatioPlotYMax = plot_features.get('RatioPlotYMax', 1.4999)  # Don't plot 1.5
        ax_ratio.set_ylim(RatioPlotYMin, RatioPlotYMax)

    # Use maximum y value from all hist datasets
    if all(isinstance(h, yoda.Scatter1D) for h in hist_data):
        max_ymax = max([h.points()[0].val(1) for h in hist_data])
    else:
        max_ymax = max([max(h.yVals()) for h in hist_data])
    if plot_features.get('YMax') is not None:
        YMax = plot_features.get('YMax')
    elif plot_features.get('LogY', 1):
        YMax = 1.7*max_ymax
    else:
        YMax = 1.1*max_ymax
    # Use minimum y value from all hist datasets
    if all(isinstance(h, yoda.Scatter1D) for h in hist_data):
        min_ymin = min([h.points()[0].val(1) for h in hist_data])
    else:
        min_ymin = min([min(h.yVals()) for h in hist_data])
    if plot_features.get('YMin') is not None:
        YMin = plot_features.get('YMin')
    elif plot_features.get('LogY', 1):
        YMin = (min_ymin/1.7 if plot_features.get('FullRange')
                else max(min_ymin/1.7, 2e-7*YMax))
    elif plot_features.get('ShowZero', 1):  # defaul ShowZero is True
        YMin = 0 if min_ymin > -1e-4 else 1.1*min_ymin
    else:
        YMin = (1.1*min_ymin if min_ymin < -1e-4 else 0 if min_ymin < 1e-4
                else 0.9*min_ymin)

    ax.set_ylim(YMin, YMax)

    # Set log scale and log tick marks frequency
    if plot_features.get('LogX'):
        ax.set_xscale('log')
        ax.xaxis.set_major_locator(mpl.ticker.LogLocator(numticks=np.inf))
    if plot_features.get('LogY', 1):
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(mpl.ticker.LogLocator(numticks=np.inf))
        ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(
            base=10.0, subs=[i for i in np.arange(0, 1, 0.1)], numticks=np.inf))

    # Set tick marks frequency given the last digit in the tick mark precision for non-log plots
    if plot_features.get('XMajorTickMarks') is not None and not plot_features.get('LogX'):
        base = plot_features.get('XMajorTickMarks')*10**(int(np.log10(XMax))-1)
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base))
    if plot_features.get('YMajorTickMarks') is not None and not plot_features.get('LogY', 1):
        base = plot_features.get('YMajorTickMarks')*10**(int(np.log10(YMax))-1)
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base))

    if plot_features.get('XMinorTickMarks') is not None and not plot_features.get('LogX'):
        ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(
            1+plot_features.get('XMinorTickMarks')))
    if plot_features.get('YMinorTickMarks') is not None and not plot_features.get('LogY', 1):
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(
            1+plot_features.get('YMinorTickMarks')))

    # Add custom ticks for x and y axes
    if plot_features.get('XCustomMajorTicks') is not None:
        ax.set_xticks(plot_features.get('XCustomMajorTicks')[::2])
        ax.set_xticklabels(plot_features.get('XCustomMajorTicks')[1::2])
        ax.set_xticks([], minor=True)  # Turn off minor xticks
        if plot_features.get('RatioPlot', 1):
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

    if plot_features.get('Legend', 1) and yoda_type == 'hist':
        handles = [AnyObject()]
        labels = [hist_features[0].get('Title', 'Data')]
        for i, _ in enumerate(hist_data[1:]):
            color = colors[i % len(colors)]
            handles.append(mpl.lines.Line2D([], [], color=color))
            labels.append(hist_features[i].get('Title', 'mc{}'.format(i+1)))
    if plot_features.get('Legend', 1) and yoda_type == 'scatter':
        handles = []
        labels = []
        for i, _ in enumerate(hist_data):
            color = colors[i % len(colors)]
            handles.append(mpl.lines.Line2D([], [], color=color))
            labels.append(hist_features[i].get('Title', 'mc{}'.format(i+1)))
    if plot_features.get('Legend', 1):
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
    if plot_features.get('RatioPlot', 1) and yoda_type == 'hist':
        fig.align_ylabels((ax, ax_ratio))
    fig.savefig(filename+'.pdf')
    fig.savefig(filename+'.png', bbox_inches="tight")
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
