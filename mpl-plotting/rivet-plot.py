"""Module creates a rivet-style plot as a pdf."""
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import data  # Import xlow, xhigh, val, errminus, errplus of data
import mc1  # Import xlow, xhigh, val, errminus, errplus of mc1

mc_error = True  # TODO: Replace this with file_style key

file_style = {'Title': 'Charged particle $p_\perp$ at 7 TeV, track $p_\perp > 500$ MeV, for $N_\mathrm{ch} \geq 1$',
              'XLabel': '$p_\perp$ [GeV]',
              'YLabel': '$1/N_\mathrm{ev} \, $1$/2\pi{}p_\perp \, \mathrm{d}\sigma/\mathrm{d}\eta\mathrm{d}p_\perp$',
              # 'ZLabel' = '...'
              # 'XLabelSep' = <distance >,
              # 'YLabelSep' = <distance >,
              # 'ZLabelSep' = <distance >,
              # 'XMajorTickMarks' = <last_digit > ,
              # 'YMajorTickMarks' = <last_digit > ,
              # 'ZMajorTickMarks' = <last_digit > ,
              # 'XMinorTickMarks' = <nticks > ,
              # 'YMinorTickMarks' = <nticks > ,
              # 'ZMinorTickMarks' = <nticks > ,
              # 'XTwosidedTicks'=<0|1>,
              # 'YTwosidedTicks'=<0|1>,
              # 'XCustomMajorTicks'=<list>
              # 'YCustomMajorTicks'=<list>
              # 'ZCustomMajorTicks'=<list>
              # 'XCustomMinorTicks'=<list>
              # 'YCustomMinorTicks'=<list>
              # 'ZCustomMinorTicks'=<list>
              # 'PlotXTickLabels'=<0|1>
              # 'RatioPlotTickLabels'=<0|1>
              'LogX': 1,
              'LogY': 1,
              # 'LogZ'=<0|1>
              'XMin': min(data.xlow),
              'XMax': 50,
              'YMin': 6e-7,
              'YMax': 3,
              # 'ZMin'=<value>
              # 'ZMax'=<value>
              # 'FullRange'=<0|1>
              # 'ShowZero'=<0|1>
              # 'NormalizeToIntegral'=<1|0>
              # 'NormalizeToSum'=<1|0>
              # 'Scale'=<factor>
              # 'Rebin'=<nbins>
              # 'PlotSize'=<xsize,ysize>
              # 'LeftMargin'=<size>
              # 'RightMargin'=<size>
              # 'TopMargin'=<size>
              # 'BottomMargin'=<size>
              # 'FrameColor'=<color>
              'Legend': 1,
              # 'CustomLegend'=<text>
              # 'LegendXPos'=<pos>
              # 'LegendYPos'=<pos>
              # 'LegendAlign'=<align>
              # 'LegendOnly'=<list>
              # 'DrawOnly'=<list>
              # 'Stack'=<list>
              # 'DrawSpecialFirst'=<0|1>
              # 'DrawFunctionFirst'=<0|1>
              # 'ConnectGaps'=<0|1>
              'RatioPlot': 1,
              # 'RatioPlotReference'=<histogram_ID>
              # 'RatioPlotMode'=<default|deviation|datamc>
              'RatioPlotYLabel': 'MC/Data',
              'RatioPlotYMin': 0.5,
              'RatioPlotYMax': 1.4999,
              # 'RatioPlotYSize'=<size>
              # 'RatioPlotErrorBandColor'=<color>
              # 'RatioPlotSameStyle'=1
              # 'MainPlot'=0
              # 'GofType'=chi2
              # 'GofReference'=<histogram_ID>
              # 'GofLegend'=<0|1>
              # 'GofFrame'=<histogram_ID>
              # 'GofFrameColor'=<colorthresholds>
              # 'ColorSeries'={rgb}{last}[rgb]{1,0.97,0.94}[rgb]{0.6,0.0,0.05}
              # 'LineStyle'=<style>
              # 'LineColor'=<color>
              # 'LineOpacity'=<opacity>
              # 'LineWidth'=<width>
              # 'LineDash'=<dashstyle>
              # 'ConnectBins'=<0|1>
              # 'ConnectGaps'=<0|1>
              # 'SmoothLine'=<0|1>
              # 'FillStyle'=<style>
              # 'FillColor'=<color>
              # 'FillOpacity'=<opacity>
              # 'HatchColor'=<color>
              'ErrorBars': 1,
              # 'ErrorBands'=<0|1>
              # 'ErrorBandColor'=<color>
              # 'ErrorBandOpacity'=<opacity>
              # 'PolyMarker'=<dotstyle>
              # 'DotSize'=<size>
              # 'DotScale'=<factor>
              # 'NormalizeToIntegral'=<1|0>
              # 'NormalizeToSum'=<1|0>
              # 'Scale'=<factor>
              # 'Rebin'=<nbins>
              # 'ErrorType'=<stat|env>
              }


def rivet_plot(yaml_file):
    # Parse yaml_file to return 1. rcparams dict, 2. plot params, 3. data dict,
    # TODO: add custom rcparams
    plt.style.use(os.path.join(sys.path[0], 'rivet.mplstyle'))

    if file_style.get('RatioPlot'):
        fig, (ax, ax_ratio) = plt.subplots(2, 1, sharex=True,
                                           gridspec_kw={'height_ratios': (2, 1)})
    else:
        fig, ax = plt.subplots(1, 1)
    ax.set_xlabel(file_style.get('XLabel'))
    ax.set_ylabel(file_style.get('YLabel'), loc='top')
    ax.set_title(file_style.get('Title'), loc='left')

    XMin = file_style.get('XMin', min(data.xlow))
    XMax = file_style.get('XMax', max(data.xhigh))
    ax.set_xlim(XMin, XMax)
    # TODO: the hist_y defaults probably need to allow more space?
    YMin = file_style.get('YMin', min(data.val))
    YMax = file_style.get('YMax', max(data.val))
    ax.set_ylim(YMin, YMax)

    # Set log scale
    if file_style.get('LogX'):
        ax.set_xscale('log')
        ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0))
    if file_style.get('LogY'):
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0))

    # hist plot
    # Data line
    ax.hlines(data.val, data.xlow, data.xhigh, 'k')
    x_data = (data.xlow + data.xhigh)/2
    y_data = data.val
    ax.plot(x_data, y_data, 'ko')
    if file_style.get('ErrorBars') and data.errminus is not None:
        ax.vlines(x_data, (data.val - data.errminus),
                  (data.val+data.errplus), 'k')
    # mc1 line
    x_mc1 = np.append(mc1.xlow, mc1.xhigh[-1])
    y_mc1 = np.insert(mc1.val, 0, mc1.val[0])
    ax.plot(x_mc1, y_mc1, 'r', drawstyle='steps-pre',
            solid_joinstyle='miter')
    if file_style.get('ErrorBars') and mc1.errminus is not None:
        ax.vlines(x_data, (mc1.val - mc1.errminus),
                  (mc1.val+mc1.errplus), 'r')

    if file_style.get('RatioPlot'):  # TODO: Check if errorbar/errorbands
        ax_ratio.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        ax_ratio.set_ylabel(file_style.get('RatioPlotYLabel', 'MC/Data'))
        RatioPlotYMin = file_style.get('RatioPlotYMin', 0.5)
        RatioPlotYMax = file_style.get('RatioPlotYMax', 1.4999)
        ax_ratio.set_ylim(RatioPlotYMin, RatioPlotYMax)

        x_ratio = x_mc1
        y_ratio = y_mc1/np.insert(y_data, 0, y_data[0])
        ax_ratio.plot(x_ratio, y_ratio, 'r', drawstyle='steps-pre', zorder=1,
                      solid_joinstyle='miter')
        ax_ratio.hlines(1, XMin, XMax, 'k', zorder=2)
        ax_ratio.plot(x_data, np.ones(len(x_data)), 'ko', zorder=3)
        ax_ratio.vlines(x_data, (mc1.val - mc1.errminus)/data.val,
                        (mc1.val + mc1.errplus)/data.val, 'r', zorder=1)

        errbandminus = np.insert((data.val - data.errminus)/data.val, 0,
                                 ((data.val - data.errminus)/data.val)[0])
        errbandplus = np.insert((data.val + data.errplus)/data.val, 0,
                                ((data.val + data.errplus)/data.val)[0])
        ax_ratio.fill_between(x_ratio, errbandminus, errbandplus,
                              step='pre', alpha=0.5, zorder=0)

    # Legend
    # TODO: Find a better way to implement the custom Data legend graphic

    if file_style.get('Legend'):
        handler_mc1 = mpl.lines.Line2D([], [], color='red')
        ax.legend([AnyObject(), handler_mc1], ['Data', 'mc1'],
                  handler_map={AnyObject: AnyObjectHandler()}, loc=[0.52, 0.77])

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
    rivet_plot(sys.argv[0])
