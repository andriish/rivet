# YAML Syntax Specification for .plot and .dat Files
The new `rivet-mkdat` command takes .plot files with a new YAML-based syntax as input and outputs the .dat files in a new format, also based on YAML. For backwards-compatibility, .plot files with the old format can also be read by `rivet-mkdat`, although there is some lacking support.

The files have sections which are formatted in the following way: 
```yaml
section name 1: 
    parameter 1: some value
    parameter 2: some other value
    ...
section name 2: 
    parameter 3: some third value
    ...
...
```
where the name of each section and its purpose is described below.

## Difference between .plot and .dat files
The structure of the .plot files and .dat files are similar, since the plot settings in a .plot file are directly added to a .dat file. There are three main differences: 
- .dat files have an additional `histogram` section, which contains the histogram data that will be plotted.
- .plot files have a `name` section which is an ID used by `rivet-mkdat` to identify which plot settings should be used for which plot. This will also be the file name of the .dat file, excluding the .dat file extension.
- There is one .plot file for each analysis (e.g., `ALICE_2010_S8624100`) but an analysis can have multiple histogram subsections (TODO not sure what these are called. Histograms?), e.g., `d11-x01-y01`. Each subsection of an analysis in a .plot file is separated by `---`.

All of these sections are described in detail below.

## plot features
The `plot features` section describes Rivet-specific plot settings. Most of these parameters already existed in the original .plot files. This section corresponds to the `BEGIN PLOT ... END PLOT` sections of the .plot files with the old format.

A list of all supported parameters are listed below.

### Camera position
```
3DElev = <angle (in degrees)>
3DAzim = <angle (in degrees)>
RatioPlot3DElev = <angle (in degrees)>
RatioPlot3DAzim = <angle (in degrees)>
```
The altitude ("latitude") and azimuthal ("longitude") angle of the "camera position" when viewing a surface plot.
If the ratio plot angles are not specified, they are set to the same value as the main plot.

TODO continue writing this list.

## style
The style is a single string which specifies which matplotlib style that should be used to change the look of a figure. The default style, called `default` (TODO: this name will be changed to e.g., `rivet-default`, since matplotlib already has a style called `default`), will create plots with the default Rivet look.

In the future, any matplotlib builin style can be specified here, e.g., `seaborn`.

### Examples
```yaml
style: default  # Apply the default rivet style. Warning: this will be renamed to e.g., rivet-default in the future
```
```yaml
style: seaborn # This does not work yet but will work in the future
```

## rcParams
This is section contains matplotlib rcParams that are passed directly to matplotlib. These settings will overwrite the corresponding setting specified by the `style` section.
A list of all available rcParams can be found in the [documentation of matplotlib](https://matplotlib.org/stable/tutorials/introductory/customizing.html#a-sample-matplotlibrc-file).

### Examples
```yaml
rcParams: 
    lines.linewidth: 2
    axes.facecolor: black
```

## histograms
This section only exists in .dat files and contains all histograms that will be plotted. Each histogram is specified by its name, which is typically the file name of the yoda file in which the histogram originally existed in. Each subsection contains a `yoda` section, which includes the histogram, printed in the yoda format. This section is often very long and not meant to be read or modified by a human. 
Additionally, the name of the histogram is specified in the `title` section and `Errorbars` specifies whether that histogram shuld plot its errorbars or not (only applied to 1D histograms).

### Examples
```yaml
histograms:
    example-yoda-file.yoda:
            ErrorBars: true
            Title: example-yoda-file
            yoda: |+
                BEGIN YODA_HISTO2D_V2 /ALEPH_2019_I1737859/backLab
                Path: /ALEPH_2019_I1737859/backLab
                ScaledBy: 3.51617440225035168e-05
                Title: ~
                Type: Histo2D
                # etc...
```

## name
This is a unique section to .plot files and specifies the path, which is a unique ID to that plot.

### Examples
```yaml
name: /ALICE_2010_S8624100/d11-x01-y01
```

# Examples
In the documentation, there is an [example .plot file](example.plot) and an [example .dat file](example.dat).