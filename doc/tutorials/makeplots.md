# make-plots

The plots produced with `rivet-mkhtml` are really rendered with the command `make-plots`, which is called under the hood. The `make-plots` command can also be used to create figures from the simple .dat text format produced by `rivet-cmphistos` directly.
The `make-plots` script is quite powerful, and includes several options to modify plotting style, adding curves or fits and more. For use with Rivet, the syntax documented in this document should be provided in the .plot file.
Internally, `make-plots` the simple text format and converts them into PostScript or PDF files by creating a LaTeX file and running `latex`, `dvips`, and maybe `ps2pdf`.

## Usage

To run `make-plots` call

```
     make-plots [options] file.dat [file2.dat ...]
```

All available options can be listed by running

```
     make-plots --help
```

### Configuration files


`make-plots` typically takes the plotting instructions and settings from the input ascii files as described in the "Input Format" chapter. It is also possible though to pass a global configuration file to `make-plots` (cf. `--help`) which allows to specify/overwrite settings for certain plots or histograms in a plot on top of what the input files specify. This could be useful if the ascii files are generated automatically (e.g. with `rivet-mkhtml` or `compare-histos`) and you still want to apply custom plotting options.

An example for this looks like:

```
    # BEGIN PLOT figures/MC_WJETS/W_mass.dat
    XMin=60.0
    XMax=100.0
    LegendXPos=0.65
    # END PLOT

    .*myLOrun.yoda/D0_2008_S7554427/d01-x01-y01::Scale=1.0

```

Here first the options in the `PLOT` section of a specific ascii file are being amended/overwritten. The second part shows how to overwrite the `Scale` property of one specific histogram line using the ID of the histogram.

## Input Format

The ascii files which can be read by `make-plots` are divided into sections. There are four types of sections which are called `PLOT`, `HISTOGRAM`, `FUNCTION`, and `SPECIAL`. Every file must contain exactly one `PLOT` section and at least one section of the other three types. There may be multiple `HISTOGRAM`, `FUNCTION`, and `SPECIAL` sections.

Empty lines and lines starting with `#` are ignored, except for the section delimiters described below.

### PLOT

The `PLOT` section starts with

```
    # BEGIN PLOT

```
and ends with


```
    # END PLOT

```

Every file must have exactly one `PLOT` section. In this section global parameters are specified, like the axis labels, the plot title, size, ??? An empty `PLOT` section is perfectly legal, though.
In this section the following parameters can be set:

#### Titles, Labels

```
    Title=<title>

```
The title of the plot.

```
    XLabel=<label>
    YLabel=<label>
    ZLabel=<label>

```
Axis labels for the x-, y-, and z-axis.


```
    XLabelSep=<distance>
    YLabelSep=<distance>
    ZLabelSep=<distance>

```
Distance between the axis label and the plot in units of `\labelsep`.

```
    XMajorTickMarks=<last_digit>
    YMajorTickMarks=<last_digit>
    ZMajorTickMarks=<last_digit>
    XMinorTickMarks=<nticks>
    YMinorTickMarks=<nticks>
    ZMinorTickMarks=<nticks>
```

`make-plots` tries to guess the distance between tickmarks automatically. If you are not satisfied with its result, you can override this by setting `<last_digit>` to 1, 2, 5, or 10, and `<nticks>` to the number of minor ticks you like. _Note_: These options are not available for logarithmic axes.

```
    XTwosidedTicks=<0|1>
    YTwosidedTicks=<0|1>

```
Draw tickmarks also on the upper and/or right side of the plot.

```
    XCustomMajorTicks=<list>
    YCustomMajorTicks=<list>
    ZCustomMajorTicks=<list>

```

To specify major ticks at arbitrary positions and/or with arbitrary labels. `<list>` is a whitespace-separated list of format `value1 <spaces_or_tabs> label1 <spaces_or_tabs> value2 <spaces_or_tabs> label2 ...`.


[//]: # TODO: allow use of YAML-style list syntax to clarify delimiters?

```
    XCustomMinorTicks=<list>
    YCustomMinorTicks=<list>
    ZCustomMinorTicks=<list>
```
To specify minor ticks at arbitrary positions. `<list>` is a tab separated list of format `value1 <tab> value2 <tab> value3 ...`.

```
    PlotXTickLabels=<0|1>
    RatioPlotTickLabels=<0|1>

```
Disable/enable plotting of the tick labels in the plot and ratio plot (useful if multiple plots are to be combined manually later).


#### Axes

```
    LogX=<0|1>
    LogY=<0|1>
    LogZ=<0|1>

```
Use a logarithmic x-, y-, or z-axis. Default is linear.

```
    XMin=<value>
    XMax=<value>
    YMin=<value>
    YMax=<value>
    ZMin=<value>
    ZMax=<value>
    FullRange=<0|1>
    ShowZero=<0|1>

```
Specify the plot range. By default the range is chosen such that all data is visible in linear plots, and the zero is visible. `ShowZero=0` suppresses plotting the zero in linear plots and thus zooms into the actual y-value range of the distribution. In logarithmic plots the automatic choice of `YMin` is limited to be not smaller than 2e-4*`YMax`, but manually you can specify any value. `FullRange=1` also overrides the 2e-4*`YMax` limit and plots the full range in y.

#### Normalization, Rebinning

```
    NormalizeToIntegral=<1|0>
    NormalizeToSum=<1|0>
    Scale=<factor>
```
Normalize all histograms to their integral, to their sum of entries, or scale them by some arbitrary factor. Normalization and scale options in the `PLOT` section override the corresponding option in the `HISTOGRAM` section. The scale factor is applied after normalization.

```
    Rebin=<nbins>
```
Rebin all histograms in this plot. Syntax and functionality is the same as for the Rebin option in the `HISTOGRAM` section.

#### Sizes and Margins

```
     PlotSize=<xsize,ysize>
```
Size in x and y direction of the plot. This can be specified in any unit LaTeX understands.

```
    LeftMargin=<size>
    RightMargin=<size>
    TopMargin=<size>
    BottomMargin=<size>
```

Distance between the plot and the paper edge.

```
    FrameColor=<color>
```
Background color for the margin around the plot.

#### Legends

```
    Legend=<0|1>
```
Display a legend in the plot.

```
    CustomLegend=<text>
```

Custom text that is added to the legend.

```
    LegendXPos=<pos>
    LegendYPos=<pos>
```
Position of the legend within the plot. Anchor point is the top left corner of the legend, so units typically range between 0.0 and 1.0.

```
    LegendAlign=<align>
```
Horizontal alignment of the legend: `LegendAlign=l` is the default and will create a left-aligned legend, while `LegendAlign=r` is right-aligned with the keys on the right hand side.

```
    LegendOnly=<list>
```
Whitespace separated list of IDs. These can be histograms or functions. The legend is only shown for the listed objects. Without this option, all plotted objects which have a title enter the legend. The legend titles are plotted in the given order, so there are cases in which it makes sense to use `LegendOnly` together with all histogram IDs. It is also possible to specify the legend order on an entry-by-entry basis using the `LegendOrder=<int>` setting for each histogram or function.


#### Plotting Options

```
    DrawOnly=<list>
```

Whitespace separated list of histogram IDs. Only the histograms in this list are plotted, even if there are more histograms defined in the file. The histograms are plotted in the given order, so there are cases in which it makes sense to use `DrawOnly` together with all histogram IDs. This is especially useful for the `Stack` option. It is also possible to specify the plotting order on a histogram-by-histogram basis using the `PlotOrder=<int>` setting for each histogram.

```
    Stack=<list>
```
Whitespace separated list of histogram IDs. The histograms will be added on top of each other. This is useful for example to compare data with background if the background has contributions from several histograms.

```
    DrawSpecialFirst=<0|1>
    DrawFunctionFirst=<0|1>
```
By default the `SPECIAL` and `FUNCTION` sections are plotted after the histograms. With these options you can override that behaviour.

```
    ConnectGaps=<0|1>
```
If error bars are disabled and you want to bridge gaps in a histogram, you can set this parameter. By default it is off. Setting it in the `PLOT` section affects all histograms, but you can also set it in the `HISTOGRAM` section for individual histograms. The local setting overrides the global setting.


#### Comparison Plots

With the

```
    RatioPlot=1
    RatioPlotReference=<histogram_ID>
```
options you can create ratio plots for two or more histograms. Note that you must specify your reference data ID. This option is used by the `compare-histos` script.

```
    RatioPlotMode=<default|deviation|datamc>
```

By default, the ratio plot displays MC/Data. You can switch to (MC-data)/uncertainty (`deviation`) or Data/MC (`datamc`) with this option.

In ratio plots the following additional options are available and work in a similar way as their regular counterparts:

```
    RatioPlotYLabel=<label>
    RatioPlotYMin=<value>
    RatioPlotYMax=<value>
    RatioPlotYSize=<size>
    RatioPlotErrorBandColor=<color>
```

By default, the reference data is plotted using a yellow error band around the central value of the ratio plot. If you would rather have it plotted in the same style as in the main plot (e.g. with black errorbars), you can specify:

```
    RatioPlotSameStyle=1
```

If you only want the ratio plot without showing the actual data distribution, you can switch off the main plot. This option implies `RatioPlot=1`:

```
    MainPlot=0
```

#### Goodness of Fit


`make-plots` can calculate the goodness of fit between histograms and display the result in the legend. It is also possible to change the color of the margin around the plot depending on the GoF. This is useful to provide a quick overview when looking at many plots.


```
    GofType=chi2
```

The type of GoF. The default is `chi2` and currently that???s the only option.


```
    GofReference=<histogram_ID>
```

specifies the reference histogram to be used for the GoF calculation. If this option is omitted, the fallback is `RatioPlotReference`.


The GoF calculation is activated by two options:

```
    GofLegend=<0|1>
    GofFrame=<histogram_ID>
```

`GofLegend` calculates the GoF for all histograms and displays the results in the legend. With `GofFrame` you can specify a single histogram for which the GoF result will be shown in the legend and used to assign a color to the plot margins. Note that `FrameColor` overrides the color choice for the margin. You can use

```
    GofFrameColor=<colorthresholds>

```

to specify the thresholds for the frame color. This option takes a list of `<threshold>:<color>` pairs, separated by whitespace. The default is `GofFrameColor=0:green 3:yellow 6:red!70`. Again, if you use `FrameColor`, this option is disabled.

#### Color Palettes for 2-dim Plots

With the option `ColorSeries` you can define a custom color palette for 2-dimensional plots. The syntax is the same as for the `\definecolorseries` command in the `xcolor` LaTeX package after the color series name, i.e. `{core-model}{method}[begin-model]{begin-spec}[end-model]{end-spec}`. For more information you can consult the [xcolor documentation](http://www.ctan.org/tex-archive/macros/latex/contrib/xcolor/xcolor.pdf). Here is an example:

```
    ColorSeries={rgb}{last}[rgb]{1,0.97,0.94}[rgb]{0.6,0.0,0.05}
```
### HISTOGRAM

The `HISTOGRAM` section starts with

```
    # BEGIN HISTOGRAM <ID>
```

and ends with

```
    # END HISTOGRAM
```

There can be more than one `HISTOGRAM` section in a file. Histograms are identified by `<ID>` which can be any string _not_ containing whitespace.

#### Data Format

Lines starting with a number (positive or negative) are interpreted as data. Each line specifies one bin. The fields in each line must be separated by tabs, not spaces (this needs to be fixes some day). For 1-dimensional histograms the format can be

```
    <lowerbinedge>  <upperbinedge>  <value>  <error>
    <lowerbinedge>  <upperbinedge>  <value>  <minuserror>  <pluserror>
```

2-dimensional histograms are supported, too. They are plotted as colormap (errors are ignored) and specified as

```
    <lowerxbinedge>  <upperxbinedge>  <lowerybinedge>  <upperybinedge>  <value>  <error>
```

#### Titles

```
    Title=<title>
```

Title of the histogram. This is used for the legend.

#### Linestyles

```
    LineStyle=<style>
```

Any linestyle that is understood by the LaTeX pstricks package, e.g. `solid`, `dotted`, `dashed`, `none`, as well as a special `dashdotted` (or `dotdashed`) linestyle which does what you might expect.

```
    LineColor=<color>
```
Color of the line. Default is black, but any color that pstricks understands can be used, including constructions like `red!70!blue!20` (for mixing colors), `{[rgb]{0.8,0,0.7}}` (for RGB-colors), `{[wave]{580}}` (for wavelengths in nm), `LineColor={[cmyk]{1,1,0,0}}` for CMYK-colors, or `[hsb]{0.5,1,1}` for HSB-colors.

```
    LineOpacity=<opacity>
```

Set the opacity of the line. Default is 1.0\. This might not work for ps output.

```
    LineWidth=<width>
```

Width of the line.

```
    LineDash=<dashstyle>
```

If `LineStyle` is set to `dashed`, you can specify the dash style with this option. Anything that is understood by pstrick???s `dash=...` option is valid. An example for a dash-dotted line is `LineDash=3pt 3pt .8pt 3pt`. You can use `LineStyle=dashdotted` or `LineStyle=dotdashed` as an abbreviation for `LineStyle=dashed` with `LineDash=3pt 3pt .8pt 3pt`.

```
    ConnectBins=<0|1>
```

Choose whether to connect adjacent bins' horizontal lines together by a vertical line on the bin edge. This is enabled by default, but you may wish to disable it when plotting reference data with error bars and point markers.

```
    ConnectGaps=<0|1>
```

If ConnectBins is enabled and you want to bridge gaps in a histogram, you can set this parameter. By default it is off. Setting it in the `PLOT` section affects all histograms, but you can also set it in the `HISTOGRAM` section for individual histograms. The local setting overrides the global setting.

```
    SmoothLine=<0|1>
```

Draw a smooth curve rather than a histogram

#### Fillstyles

```
    FillStyle=<style>
    FillColor=<color>
```

To fill the area below a histogram, set `FillStyle` and `FillColor` to something pstricks understands. Examples for the style are `solid` or `vlines`. See `LineColor` for examples of color definitions.

```
    FillOpacity=<opacity>
```

Set the opacity of the solid fillcolor. Default is 1.0\. This might not work for ps output.

```
    HatchColor=<color>
```

The color of a hatch pattern used for filling the area below a histogram. This is used for example when you use `vlines` as style.

#### Data Points

```
    ErrorBars=<0|1>
```

Turn on error bars.

```
    ErrorBands=<0|1>
    ErrorBandColor=<color>
```
Turn on error bands and set their color (see `LineColor` for a description of color definitions).

```
    ErrorBandOpacity=<opacity>
```

Set the opacity of the error band. Default is 1.0\. This might not work for ps output.

```
    PolyMarker=<dotstyle>
```

The marker style of the points. Any dot style which is understood by pstricks is valid, e.g. `*`, `o`, `triangle`, `diamond`, ???

```
    DotSize=<size>
    DotScale=<factor>
```

The size of the markers. With `DotSize` you can specify the absolute size, e.g. in units of `pt`, while `DotScale` is a relative measure with respect to the default size.

#### Normalization, Rebinning

```
    NormalizeToIntegral=<1|0>
    NormalizeToSum=<1|0>
    Scale=<factor>
```

Normalize the histogram to the integral, to the sum of entries, or scale it by some arbitrary factor. If normalization and a scale factor are given, the scale factor is applied after normalization. This is useful for stacking histograms when the ratios are known.

```
    Rebin=<nbins>
    ErrorType=<stat|env>
```

Rebin the histogram. Starting with the lowest bin `<nbins>` bins are combined into a new bin. If the number of bins in the histogram is not a multiple of `<nbins>`, the remaining bins at the upper histogram end are silently ignored (i.e. if the original histogram has 10 bins and `<nbins>` is 3, the plotted histogram shows three bins combining the bins 1???9 of the original histogram). The treatment of the errors is determined by the given ErrorType: `stat` (default) assumes the errors are of statistical nature and combines them in quadrature sum, while `env` allows to treat errors as envelope of various uncertainty runs which are combined linearly.

### FUNCTION

`make-plots` can draw arbitrary functions. These functions are defined as python code sniplets which are evaluated by `make-plots`. The code sniplet must come after all other options in a `FUNCTION` section and are preceded by `Code=` on a single line. An example `FUNCTION` section might look like this:

```
    # BEGIN FUNCTION f_cc
    LineColor=red
    Code=
    p0=16.4
    p1=1.25
    p2=0.9832
    from scipy.special import erf
    x-=0.5
    if x<=0:
        return 0
    else:
        return .5*p2*(1.+erf( (x-p0)/sqrt(x*p1) ))
    # END FUNCTION
```

#### Common Options with HISTOGRAM

The following options have the same meaning as in the `HISTOGRAM` section:

```
    Title=<title>
    LineStyle=<style>
    LineColor=<color>
    LineWidth=<width>
    LineDash=<dashstyle>
    FillStyle=<style>
    FillColor=<color>
    HatchColor=<color>
```

#### Function Range

You can limit the plot range of functions by specifying

```
    XMin=<value>
    XMax=<value>
```

### SPECIAL

The `SPECIAL` sections are used to include any custom pstricks code. This is useful for drawing arrows and lines, put text at any position into the plot, etc. The default coordinate system is defined to be `(0,0)` at the lower left and `(1,1)` at the upper right corner of the plot. By putting the `\physicscoor` command in front of a coordinate pair, these coordinates are interpreted not in the pstricks coordinate system, but in the physics coordinate system of the plot, which is useful e.g. for marking cut values in a plot. Similar `\physicsxcoor` and `\physicsycoor` commands exist which will only treat the x or y coordinate respectively as being in physics units.

*Hint:* If you want to clip your `SPECIAL` code to the plot area, you can use

```
    \psclip{\psframe[linewidth=0, linestyle=none](0,0)(1,1)}
       ...
    \endpsclip
```

An example of a `SPECIAL` section might look like this:

```
    # BEGIN SPECIAL
    \psclip{\psframe[linewidth=0, linestyle=none](0,0)(1,1)}
    \psline[linewidth=1.2pt,linecolor=red]{<-}\physicscoor(2.83,2)\physicscoor(2.83,18)
    \uput{4pt}[180]{0}\physicscoor(2.83,12){observed}
    \psline[linewidth=0.8pt,linecolor=red,linestyle=dashed]\physicscoor( 3.17,0)\physicscoor( 3.17,28.14)
    \psline[linewidth=0.8pt,linecolor=red,linestyle=dashed]\physicscoor(-3.59,0)\physicscoor(-3.59,28.14)
    \endpsclip
    # END SPECIAL
```
