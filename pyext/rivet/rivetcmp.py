#! /usr/bin/env python

"""\
Library for extracting data from YODA files and generate histogram comparison plots.

ENVIRONMENT:
 * RIVET_ANALYSIS_PATH: list of paths to be searched for plugin
     analysis libraries at runtime
 * RIVET_DATA_PATH: list of paths to be searched for data files
"""

from __future__ import print_function
import rivet, yoda, sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize, integrate
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import logging

rivet.util.check_python_version()
rivet.util.set_process_name(os.path.basename(__file__))


class Plot(dict):
    "A tiny Plot object to help keeping track of PLOT file contents. Option to write out"
    def __repr__(self):
        return "# BEGIN PLOT\n" + "\n".join("%s=%s" % (k,v) for k,v in self.items()) + "\n# END PLOT\n\n"


def sanitiseString(s):
    #s = s.replace('_','\\_')
    #s = s.replace('^','\\^{}')
    #s = s.replace('$','\\$')
    s = s.replace('#','\\#')
    s = s.replace('%','\\%')
    return s

def getHistos(args,filelist):
    """Loop over all input files. Only use the first occurrence of any REF-histogram
    and the first occurrence in each MC file for every MC-histogram."""
    refhistos, mchistos = {}, {}
    for infile in filelist:
        mchistos.setdefault(infile, {})
        analysisobjects = yoda.read(infile, patterns=args.PATHPATTERNS, unpatterns=args.PATHUNPATTERNS)
        #print(analysisobjects)
        for path, ao in analysisobjects.items():
            ## We can't plot non-histograms yet
            # TODO: support counter plotting with a faked x (or y) position and forced plot width/height
            if ao.type() not in ("Counter",
                               "Histo1D",
                               "Histo2D",
                               "Profile1D",
                               "Profile2D",
                               "Scatter1D",
                               "Scatter2D",
                               "Scatter3D"):
                continue

            ## Make a path object and ensure the path is in standard form
            try:
                aop = rivet.AOPath(path)
            except Exception as e:
                #print(e)
                print("Found analysis object with non-standard path structure:", path, "... skipping")
                continue

            ## We don't plot data objects with path components hidden by an underscore prefix
            if aop.istmp() or aop.israw():
                continue

            ## Add it to the ref or mc paths, if this path isn't already known
            basepath = aop.basepath(keepref=False)
            if aop.isref() and basepath not in refhistos:
                ao.path = aop.varpath(keepref=False, defaultvarid="0")
                refhistos[basepath] = ao
            else: #if basepath not in mchistos[infile]:
                mchistos[infile].setdefault(basepath, {})[aop.varid("0")] = ao

    return refhistos, mchistos


def getHistoLists(mchistos):
    """Make a list of all histograms and 2d histos to be plotted,
    as well as a unique list of all analyses"""
    hpaths, h2ds = [], []
    for aos in mchistos.values():
        for p in aos.keys():
            ps = rivet.stripOptions(p)
            if ps and ps not in hpaths:
                hpaths.append(ps)
            firstaop = aos[p][sorted(aos[p].keys())[0]]
            # TODO: Would be nicer to test via isHisto and dim or similar, or yoda.Scatter/Histo/Profile base classes
            if type(firstaop) in (yoda.Histo2D, yoda.Profile2D) and ps not in h2ds:
                h2ds.append(ps)

    # Unique list of analyses
    anas = list(set([x.split("/")[1] for x in hpaths]))
    return hpaths, h2ds, anas

def getRivetRefData(args,anas=None):
    "Find all Rivet reference data files"
    refhistos = {}
    rivet_data_dirs = rivet.getAnalysisRefPaths()
    dirlist = []
    for d in rivet_data_dirs:
        if anas is None:
            import glob
            dirlist.append(glob.glob(os.path.join(d, '*.yoda')))
        else:
            dirlist.append([os.path.join(d, a+'.yoda') for a in anas])
    for filelist in dirlist:
        # TODO: delegate to getHistos?
        for infile in filelist:
            analysisobjects = yoda.read(infile, patterns=args.PATHPATTERNS, unpatterns=args.PATHUNPATTERNS)
            for path, ao in analysisobjects.items():
                aop = rivet.AOPath(ao.path())
                if aop.isref():
                    ao.setPath(aop.basepath(keepref=False))
                    refhistos[ao.path()] = ao
    return refhistos

def parseArgs(args):
    """Look at the argument list and split it at colons, in order to separate
    the file names from the plotting options. Store the file names and
    file specific plotting options."""
    filelist = []
    plotoptions = {}
    for a in args:
        asplit = a.split(':')
        path = asplit[0]
        filelist.append(path)
        plotoptions[path] = []
        has_title = False
        for i in range(1, len(asplit)):
            ## Add 'Title' if there is no = sign before math mode
            if '=' not in asplit[i] or ('$' in asplit[i] and asplit[i].index('$') < asplit[i].index('=')):
                asplit[i] = 'Title=%s' % asplit[i]
            if asplit[i].startswith('Title='):
                has_title = True
            plotoptions[path].append(asplit[i])
        if path != "PLOT" and not has_title:
            plotoptions[path].append('Title=%s' % sanitiseString(os.path.basename( os.path.splitext(path)[0] )) )
    return filelist, plotoptions


def setOptions(ao, options):
    "Set arbitrary annotations"
    for opt in options:
        key, val = opt.split('=', 1)
        ao.setAnnotation(key, val)


# TODO: move to rivet.utils
def mkoutdir(outdir):
    "Function to make output directories"
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except:
            msg = "Can't make output directory '%s'" % outdir
            raise Exception(msg)
    if not os.access(outdir, os.W_OK):
        msg = "Can't write to output directory '%s'" % outdir
        raise Exception(msg)


def mkOutput(hpath, aos, plot=None, special=None):
    """
    Make the .dat file string. We can't use "yoda.writeFLAT(anaobjects, 'foobar.dat')"
    because the PLOT and SPECIAL blocks don't have a corresponding analysis object.
    """
    output = ''

    if plot is not None:
        output += str(plot)

    if special is not None:
        output += "\n"
        output += "# BEGIN SPECIAL %s\n" % hpath
        output += special
        output += "# END SPECIAL\n\n"

    from io import StringIO
    sio = StringIO()
    yoda.writeFLAT(aos, sio)
    output += sio.getvalue()

    return output


# Construct a plot object and possibly a special object for a figure
def plotSpec(hpath, is2d, plotparser, plotoptions, args=None):
    ## A Plot object to represent the PLOT file contents
    plot = Plot()
    plot['Showratio'] = '1'
    if args is not None:
        showratio = args.RATIO and not is2d
    else:
        showratio = not is2d
    if not showratio:
        plot['Showratio'] = '0'
    if not is2d:
        plot['Legend'] = '1'
        plot['LogY'] = '1'
    headers = plotparser.getHeaders(hpath)
    if headers:
        plot.update(headers)
    if "PLOT" in plotoptions:
        for key_val in plotoptions["PLOT"]:
            key, val = [s.strip() for s in key_val.split("=", 1)]
            if 'ReplaceOption' in key_val:
                opt_group = key_val.split("=", 2)
                key, val = "=".join(opt_group[:2]), opt_group[2]
            plot[key] = val
    if args is not None:
        if args.REMOVE_OPTIONS:
            plot['RemoveOptions'] = '1'
        if args.LINEAR:
            plot['LogY'] = '0'
        if args.NOPLOTTITLE:
            plot['Title'] = ''
        if showratio and args.RATIO_DEVIATION:
            plot['RatioPlotMode'] = 'deviation'
        if args.STYLE == 'talk':
            plot['PlotSize'] = '8,6'
        elif args.STYLE == 'bw' and showratio:
            plot['RatioPlotErrorBandColor'] = 'black!10'


    ## Get a special object, if there is one for this path
    special = plotparser.getSpecial(hpath)

    return plot, special

# Extract coordinates and errors from a yoda scatter.
def getXY(h_mc):
    h = yoda.mkScatter(h_mc)
    x = [b.x() for b in h]
    y = [b.y() for b in h]
    xem = [b.xErrs()[0] for b in h]
    xep = [b.xErrs()[1] for b in h]
    yem = [b.yErrs()[0] for b in h]
    yep = [b.yErrs()[1] for b in h]
    return x, y, [xem, xep], [yem,yep]

def getXYZ(h_mc):
    h = yoda.mkScatter(h_mc)
    x = [b.x() for b in h]
    y = [b.y() for b in h]
    z = [b.z() for b in h]
    xem = [b.xErrs()[0] for b in h]
    xep = [b.xErrs()[1] for b in h]
    yem = [b.yErrs()[0] for b in h]
    yep = [b.yErrs()[1] for b in h]
    zem = [b.zErrs()[0] for b in h]
    zep = [b.zErrs()[1] for b in h]
    return x, y, z, [xem, xep], [yem,yep], [zem, zep]

# Make a matplotlib canvas with or without room for ratio.
def makefig(ratio=True):
    fig = plt.figure()
    if ratio:
        spec = gridspec.GridSpec(ncols=1,nrows=2,height_ratios=[3,1])
        ax1 = fig.add_subplot(spec[0])
        ax2 = fig.add_subplot(spec[1],sharex=ax1)
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.subplots_adjust(wspace=0, hspace=0)
        return fig, ax1, ax2
    else:
        ax1 = fig.add_subplot(111)
        return fig, ax1, []

## Add a scatter to a matplotlib AxesSubplot object.
def addToAx(ax, data, plotoptions, colIter, errorbar=True, isData=False, forceTitle=''):
    # Extract file specific plot options
    lineStyle = '-'
    title = ''
    if isData:
        lineStyle = 'o'
        title = 'Data'
        lineColor = 'black'
        # Force data color to white if background is black.
        if plt.rcParams['axes.facecolor'] == 'black':
            lineColor = 'white'
    else:
        lineColor = next(colIter)
        for opt in plotoptions:
            if 'Title=' in opt:
                title = opt.replace('Title=','')
            if 'LineStyle=' in opt:
                lineStyle = opt.replace('LineStyle=','')
            if 'LineColor=' in opt:
                lineColor = opt.replace('LineColor=','')
    x, y, xe, ye = getXY(data)
    if forceTitle is not '':
        title = forceTitle
    if isData:
        ax.errorbar(x, y, xerr=xe, yerr=ye, fmt=lineStyle, color=lineColor,label=title)
    else:
        if errorbar:
            ax.errorbar(x, y, yerr=ye, drawstyle='steps-mid',fmt=lineStyle, color=lineColor, label=title)
        else:
            ax.errorbar(x, y, drawstyle='steps-mid',fmt=lineStyle, color=lineColor, label=title)
    return colIter

## Add a 3d scatter to a matplotlib AxesSubplot object.
def add3dToAx(hpath, ax, data, plot, plotoptions, cmap):
    # Extract file specific plot options
    title = ''
    plotType = 'Trisurf'
    for opt in plotoptions:
        if 'Title=' in opt:
            title = opt.replace('Title=','')
    for p in plot:
        opt = plot[p]
        opt = opt.replace("text","mathrm")
        opt = opt.replace("~"," ")
        if p == '3DPlotType':
            plotType = opt.replace('3DPlotType=','')
    x, y, z, xe, ye, ze = getXYZ(data)
    if plotType == 'Scatter':
        return ax.scatter3D(x, y, z, c = z, cmap=cmap)
    return ax.plot_trisurf(x, y, z, cmap=cmap)

## Calculate the error on the ratio.
def ratioError(ratio, yde, yd, yne, yn, ratioMode):
    if ratioMode == "default":
        rep = [rr * np.sqrt((y1e/y1)**2 + (y2e/y2)**2) for rr, y1e, y1, y2e, y2 in zip(ratio,yde[0],yd,yne[0],yn)]
        rem = [rr * np.sqrt((y1e/y1)**2 + (y2e/y2)**2) for rr, y1e, y1, y2e, y2 in zip(ratio,yde[1],yd,yne[1],yn)]
        return [rep, rem]
    elif ratioMode == "deviation": # TODO: Fix this and other options.
        rep = [rr * np.sqrt((y1e[0]/y1)**2 + (y2e[0]/y2)**2) for rr, y1e, y1, y2e, y2 in zip(ratio,yde,yd,yne,yn)]
        rem = [rr * np.sqrt((y1e[1]/y1)**2 + (y2e[1]/y2)**2) for rr, y1e, y1, y2e, y2 in zip(ratio,yde,yd,yne,yn)]
        return [rep, rem]

## Add a ratio to a matplotlib AxesSubplot object.
def addRatioToAx(hpath, ax, num, denom, plotoptions, colIter, errorbar=True, isData=False):
    if isData:
        lineStyle = 'o'
        title = 'Data'
        lineColor = 'black'
        # Force data color to white if background is black.
        if plt.rcParams['axes.facecolor'] == 'black':
            lineColor = 'white'
        ratioMode = 'default'

    else: 
        # Extract file specific plot options
        lineStyle = '-'
        lineColor = next(colIter)
        ratioMode = 'default'
        for opt in plotoptions:
            if 'LineStyle=' in opt:
                lineStyle = opt.replace('LineStyle=','')
            if 'LineColor=' in opt:
                lineColor = opt.replace('LineColor=','')
            if 'RatioPlotMode' in opt:
                ratioMode = opt.replace('RatioPlotMode=','')
                modes = ['default','deviation','datamc']
                if ratioMode not in modes:
                    print('Ratio error in '+hpath+'. Reverting to default. Choose from:')
                    print(modes)
                    ratioMode = 'default'
    xn, yn, xne, yne = getXY(num)
    xd, yd, xde, yde = getXY(denom)
    ratio = []
    re = []
    if ratioMode == 'default':
        ratio = [y1/y2 for y1, y2 in zip(yn,yd)]
    elif ratioMode == 'deviation':
        ratio = [(y1 - y2)/y2 for y1, y2 in zip(yn,yd)]
    if errorbar:
        re = ratioError(ratio, yde, yd, yne, yn, ratioMode)
    try:
        if isData:
            ax.errorbar(xd, ratio, xerr=xde, yerr=re, fmt=lineStyle, color=lineColor,label=title)
        else:
            if errorbar:
                ax.errorbar(xd, ratio, yerr=re, fmt=lineStyle, color=lineColor, drawStyle='steps-mid')
            else:
                ax.plot(xd, ratio, lineStyle=lineStyle, color=lineColor)
    except Exception as e:
        print("Ratio error in "+hpath)
        print(e)
    return colIter

# Add attributes read from the .plot file to the ax objects.
def addPlotAttributes(hpath, plot, is2d, ax1, ax2 = []):
    # Use the most normal .PLOT options
    # Default is logarithmic y-axis.
    if not is2d:
        ax1.set_yscale('log')
    # Catch exception because .PLOT files are error prone.
    showratio = True
    if plot['Showratio'] == '0':
        showratio = False
    try:
        for p in plot:
            opt = plot[p]
            opt = opt.replace("text","mathrm")
            opt = opt.replace("~"," ")
            if p == 'LogX' and int(opt.replace('LogX=','')) == 1:
                ax1.set_xscale('log')
            if p == 'LogY' and int(opt.replace('LogY=','')) == 0:
                ax1.set_yscale('linear')
            if p == 'LogZ' and is2d and int(opt.replace('LogZ=','')) == 1:
                ax1.set_zscale('log')
            if p == 'Title':
                ax1.set_title(opt)
            if p == 'XLabel':
                if showratio:
                    ax2.set_xlabel(opt,x=0.97)
                else:
                    ax1.set_xlabel(opt,x=0.97)
            if p == 'YLabel':
                ax1.set_ylabel(opt,y=0.97)
            if p == 'ZLabel' and is2d:
                ax1.set_zlabel(opt,y=0.97)
            if p == 'Legend' and opt.replace('Legend=','') == True:
                ax1.legend()
            if p == 'XMin':
                xMin = float(opt.replace('XMin=',''))
                if is2d:
                    ax1.set_xlim3d(xmin=xMin)
                else:
                    ax1.set_xlim(xmin=xMin)
                if showratio:
                    ax2.set_xlim(xmin=xMin)
            if p == 'XMax':
                xMax = float(opt.replace('XMax=',''))
                if is2d:
                    ax1.set_xlim3d(xmax=xMax)
                else:
                    ax1.set_xlim(xmax=xMax)
                if showratio:
                    ax2.set_xlim(xmax=xMax)
            if p == 'YMin':
                yMin = float(opt.replace('YMin=',''))
                ax1.set_ylim(ymin=yMin)
            if p == 'YMax':
                yMax = float(opt.replace('YMax=',''))
                ax1.set_ylim(ymax=yMax)
            if p == 'Legend' and bool(opt.replace('Legend=','')) == True:
                # Override matplotlibs silly internal legend ordering.
                handles, labels = ax1.get_legend_handles_labels()
                labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
                ax1.legend(handles, labels)
    except Exception as e:
        logging.error("Plotting error, check your .plot file. Plot: "+hpath)
        logging.error(str(e))
