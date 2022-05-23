import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yoda
import yoda_plot
import math
import re

def plot_1Dhist(hist_data, hist_features, yaml_dicts, filename):
    """Plot the Scatter1D and Scatter2D objects using Rivet styles.

    Parameters
    ----------
    hist_data : list[yoda.Scatter2D]
        All histograms that will be plotted, where the reference hist is the first element.
    hist_features : list[dict]
        Plot settings for each histogram.
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
    
    # TODO there should be a cleaner way, but \textrm is generally
    # not working in math mode. Naively check if there is a regex match
    # to a \textrm command in math mode to catch some exceptions... 
    if 'Title' in plot_features.keys() and re.search("^.*\$.*textrm.*\$.*$", plot_features['Title']):
      plot_features['Title'] = plot_features['Title'].replace("\\textrm", "\\mathrm")
    if 'XLabel' in plot_features.keys() and re.search("^.*\$.*textrm.*\$.*$", plot_features['XLabel']):
      plot_features['XLabel'] = plot_features['XLabel'].replace("\\textrm", "\\mathrm")
    if 'YLabel' in plot_features.keys() and re.search("^.*\$.*textrm.*\$.*$", plot_features['YLabel']):
      plot_features['YLabel'] = plot_features['YLabel'].replace("\\textrm", "\\mathrm")

    yoda_type = 'hist' if isinstance(hist_data[0], yoda.Scatter2D) else 'scatter'
    
    ax_format = {}  # Stores the items for formatting the axes in a dict

    # Create fig and axes
    if plot_features.get('RatioPlot', 1) and yoda_type == 'hist':
        fig, (ax, ax_ratio) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': (2, 1)})
    else:
        fig, ax = plt.subplots(1, 1)

    plt.rcParams['xtick.top'] = plot_features.get('XTwosidedTicks', True)
    plt.rcParams['ytick.right'] = plot_features.get('YTwosidedTicks', True)

    # Set plot lims
    if yoda_type == 'hist':
        XMin = plot_features.get('XMin', min([h.xMin() for h in hist_data]))
        XMax = plot_features.get('XMax', max([h.xMax() for h in hist_data]))
        ax_format['xlim'] = (XMin, XMax)

    if plot_features.get('RatioPlot', 1) and yoda_type == 'hist':
        # Ratio plot has y range of 0.5 to 1.5
        ax_ratio.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        ax_ratio.set_ylabel(plot_features.get('RatioPlotYLabel', 'MC/Data'))
        RatioPlotYMin = plot_features.get('RatioPlotYMin', 0.5)
        RatioPlotYMax = plot_features.get('RatioPlotYMax', 1.4999)  # Don't plot 1.5
        ax_ratio.set_ylim(RatioPlotYMin, RatioPlotYMax)

    # set maximum Y value from all hist datasets 
    if yoda_type == 'scatter':
        max_ymax = max([h.points()[0].val(1) for h in hist_data])
    else:
        max_ymax = max([max(h.yVals()) for h in hist_data])
    if plot_features.get('YMax') is not None:
        YMax = plot_features.get('YMax')
    elif plot_features.get('LogY', 1):
        # round off highest number in the histograms to next power of 10
        YMax = 10**(math.ceil(math.log10(max_ymax)))
    else:
        YMax = 1.1*max_ymax
    
    # Use minimum y value from all hist datasets
    if yoda_type == 'scatter':
        min_ymin = min([h.points()[0].val(1) for h in hist_data]) # TO DO -- where does this come from??
    else:
        min_ymin = min([min(h.yVals()) for h in hist_data])
    if plot_features.get('YMin') is not None:
        YMin = plot_features.get('YMin')
    elif plot_features.get('LogY', 1):
        # round off lowest number in the histograms to lower power of 10
        YMin = 10**(math.floor(math.log10(min_ymin))) if min_ymin !=0 else 2e-7*YMax # TODO: come up with a better solution to deal with min_ymin=0
    elif plot_features.get('ShowZero', 1):  # default ShowZero is True
        YMin = 0 if min_ymin > -1e-4 else 1.1*min_ymin
    else:
        YMin = (1.1*min_ymin if min_ymin < -1e-4 else 0 if min_ymin < 1e-4
                else 0.9*min_ymin)

    ax_format['ylim'] = (YMin, YMax)
    ax_format['logx'] = plot_features.get('LogX')
    ax_format['logy'] = plot_features.get('LogY', 1)
    if ax_format['logy']:
        ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(
            base=10.0, subs=[i for i in np.arange(0, 1, 0.1)], numticks=np.inf))
    ax_format['xmajor_ticks'] = plot_features.get('XMajorTickMarks')
    ax_format['ymajor_ticks'] = plot_features.get('YMajorTickMarks')
    ax_format['xminor_ticks'] = plot_features.get('XMinorTickMarks')
    ax_format['yminor_ticks'] = plot_features.get('YMinorTickMarks')
    ax_format['xcustom_major_ticks'] = plot_features.get('XMajorTickMarks')
    #if plot_features.get('nRatioTicks'): print("NRATIOTICKS: ", plot_features.get('nRatioTicks'))

    if plot_features.get('XCustomMajorTicks') is not None:
        ax_format['xcustom_major_ticks'] = plot_features.get('XMajorTickMarks')
        if plot_features.get('RatioPlot', 1):
            ax_ratio.set_xticks([], minor=True)

            
    ax_format['ycustom_major_ticks'] = plot_features.get('YMajorTickMarks')
    ax_format['ycustom_minor_ticks'] = plot_features.get('YMajorTickMarks')
    ax_format['ycustom_minor_ticks'] = plot_features.get('YMajorTickMarks')
    ax_format['plot_xticklabels'] = plot_features.get('PlotXTickLabels') 
        
    plot_errorbars = [h.get('ErrorBars', 1) for h in hist_features]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # Plot histogram and ratio using Yoda function
    if yoda_type == 'scatter':
        yoda_plot.plot_scatter1D(hist_data, ax, colors, **ax_format)
    else:
        yoda_plot.plot_hist(hist_data, True, ax, plot_errorbars, colors, **ax_format)
        if plot_features.get('RatioPlot', 1):
            yoda_plot.plot_ratio(hist_data, ax_ratio, plot_errorbars, plot_features.get('ErrorBands'), colors)

    # Create legend
    if plot_features.get('Legend', 1) and yoda_type == 'hist':
        handles = [AnyObject()]
        labels = [hist_features[0].get('Title', 'Data')]
        for i, _ in enumerate(hist_data[1:]):
            color = colors[i % len(colors)]
            handles.append(mpl.lines.Line2D([], [], color=color))
            labels.append(hist_features[i+1].get('Title', 'mc{}'.format(i+1)).replace(".yoda",""))
    if plot_features.get('Legend', 1) and yoda_type == 'scatter':
        handles = []
        labels = []
        for i, _ in enumerate(hist_data):
            color = colors[i % len(colors)]
            handles.append(mpl.lines.Line2D([], [], color=color))
            labels.append(hist_features[i].get('Title', 'mc{}'.format(i+1)).replace(".yoda",""))
      
    
  
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

    # Set text labels on axes
    if plot_features.get('RatioPlot', 1) and yoda_type == 'hist':
        ax_ratio.set_xlabel(plot_features.get('XLabel'))
    else:
        ax.set_xlabel(plot_features.get('XLabel'))
    ax.set_ylabel(plot_features.get('YLabel'), loc='top')
    ax.set_title(plot_features.get('Title'), loc='left')
      
    if plot_features.get('RatioPlot', 1) and yoda_type == 'hist':
        fig.align_ylabels((ax, ax_ratio))

    # set number of minor ticks between major ticks in ratio plot
    if plot_features.get('RatioPlot') == 1:
        nMinorTicks = int(plot_features.get('nRatioTicks')) if plot_features.get('nRatioTicks') else 1
        if nMinorTicks == 0:
            divider = 1
        else:
            divider = 0.1/(nMinorTicks+1)
        ax_ratio.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(divider))

    # remove trailing decimal zeroes for ticks, e.g. display "20" instead of "20.0"
    if not plot_features.get('LogY', 1): plt.gca().xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%g'))
    # maximum number of x-axis labels       
    plt.gca().xaxis.set_major_locator(mpl.ticker.AutoLocator())

    return fig,ax

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
