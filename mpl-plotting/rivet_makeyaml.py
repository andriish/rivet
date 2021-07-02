from __future__ import print_function
import rivet, yoda
import os, glob, io
import yamlparser 
# TODO: move all yaml-related things to a yamlio.py file so that if parser needs to be changed, it is only in that file?
from ruamel import yaml

# TODO: add more descriptive docstrings to all functions.

# TODO: remove these variables once the names have been properly decided
histogram_str_name = 'flat'


# This class, function and call to add_representer makes it so that all objects with type literal will be printed to a .yaml file as a string block. 
class literal(str):
    """A small wrapper class used to print histograms as multiline strings in a .yaml file."""
    pass

def literal_presenter(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(literal, literal_presenter)


def sanitiseString(s):
    s = s.replace('#','\\#')
    s = s.replace('%','\\%')
    return s


def _parse_args(args):
    """Look at the argument list and split it at colons, in order to separate
    the file names from the plotting options. Store the file names and
    file specific plotting options.
    
    Parameters
    ----------
    args : list[str]
        List of arguments which were previously passed to `rivet-cmphistos`. Will be ['filename.yoda:key=value', ..., 'PLOT:key=value:key=value'
    
    Returns
    -------
    filelist : list[str]
        Raw names of the files, i.e. the first part of each string in args.
    filenames : list[str]
        Names of the files. If Name=value is passed as a plot option after a file name, this will become the filename. Otherwise, it will use the same value as in filelist.
    plotoptions : dict[str, list[str]]
        Dictionary of plot options. 
        The key will be the file name (i.e., same value as in filenames) and the value will be a list of strings with plot options. 
        One of the keys will also be PLOT (if it was passed in as an argument to args, which contains all plot options that will be applied to all input files.
    
    Examples
    --------
    >>> _parse_args(['mc1.yoda:Title=example title:Name=example name 1', 'mc2.yoda', 'PLOT:LogX=1'])
    (['mc1.yoda', 'mc2.yoda'],
     ['example name 1', 'mc2.yoda'],
     {
         'example name 1': ['Title=example title', 'Name=example name 1'],
         'mc2.yoda': ['Title=mc2'],
         'PLOT': ['LogX=1']
    })
    
    Note
    ----
    Some matplotlib line styles contain ':', which would not work with current code. TODO:  change delimiter?
    """
    # TODO: remove filenames since they exist as keys in plotoptions?
    filelist = []
    filenames = []
    plotoptions = {}
    for a in args:
        asplit = a.split(':')
        path = asplit[0]
        if path != "PLOT":
            filelist.append(path)
            filenames.append(path)
        plotoptions[path] = {}
        has_title = False
        has_name = ""
        for i in range(1, len(asplit)):
            ## Add 'Title' if there is no = sign before math mode
            if '=' not in asplit[i] or ('$' in asplit[i] and asplit[i].index('$') < asplit[i].index('=')):
                asplit[i] = 'Title=%s' % asplit[i]
            if asplit[i].startswith('Title='):
                has_title = True
            key, value = asplit[i].split('=', 1)
            plotoptions[path][key] = value
            if asplit[i].startswith('Name=') and path != "PLOT":
                has_name = asplit[i].split('=', 1)[1]
                filenames[-1] = has_name
        if has_name != "":
            plotoptions[has_name] = plotoptions[path]
            del plotoptions[path]
        if path != "PLOT" and not has_title:
            plotoptions[has_name if has_name != "" else path]['Title'] = sanitiseString(os.path.basename( os.path.splitext(path)[0] ))
    return filelist, filenames, plotoptions


def _get_histos(filelist, filenames, plotoptions, path_patterns, path_unpatterns):
    """Loop over all input files. Only use the first occurrence of any REF-histogram
    and the first occurrence in each MC file for every MC-histogram."""
    refhistos, mchistos = {}, {}
    for infile, inname in zip(filelist, filenames):
        mchistos.setdefault(inname, {})
        analysisobjects = yoda.read(infile, patterns=path_patterns, unpatterns=path_unpatterns)
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
                print("Found analysis object with non-standard path structure:", path, "... skipping")
                continue

            ## We don't plot data objects with path components hidden by an underscore prefix
            if aop.istmp() or aop.israw():
                continue

            ## Add it to the ref or mc paths, if this path isn't already known
            basepath = aop.basepath(keepref=False)
            defaultWeightName = plotoptions[inname].get('DefaultWeight', '0')
            if aop.isref() and basepath not in refhistos:
                ao.setPath(aop.varpath(keepref=False, defaultvarid=defaultWeightName))
                refhistos[basepath] = ao
            else: #if basepath not in mchistos[infile]:
                mchistos[inname].setdefault(basepath, {})[aop.varid(defaultWeightName)] = ao

    return refhistos, mchistos


def _get_rivet_ref_data(anas, path_patterns, path_unpatterns):
    """Find all Rivet reference data files"""
    refhistos = {}
    rivet_data_dirs = rivet.getAnalysisRefPaths()
    dirlist = []
    for d in rivet_data_dirs:
        if anas is None:
            dirlist.append(glob.glob(os.path.join(d, '*.yoda*')))
        else:
            #dirlist.append([os.path.join(d, a+'.yoda*') for a in anas])
            for a in anas:
                res = glob.glob(os.path.join(d, a+'.yoda'))
                if len(res) == 0:
                    res = glob.glob(os.path.join(d, a+'.yoda.gz'))
                if len(res) != 0:
                    dirlist.append(res)
    for filelist in dirlist:
        # TODO: delegate to _get_histos?
        for infile in filelist:
            analysisobjects = yoda.read(infile, patterns=path_patterns, unpatterns=path_unpatterns)
            for path, ao in analysisobjects.items():
                aop = rivet.AOPath(ao.path())
                if aop.isref():
                    ao.setPath(aop.basepath(keepref=False))
                    refhistos[ao.path()] = ao
    return refhistos


def _make_output(plot_id, plotdirs, config_files, mchistos, refhistos, reftitle, filelist, filenames, plotoptions):
    """
    Create output dictionary for the plot_id.
    
    Parameters
    ----------
    plot_id : str
        ID, usually of the format AnalysisID/HistogramID.
    plotdirs : list[str]
        All directories to look for .plot files at.
    config_files : list[str]
        Additional plot settings that will be applied to all figures.
    mchistos : dict
        Dictionary of the Monte Carlo YODA histograms.
        The structure is {filename: {plot_id: {"0": yoda_histogram1, "1": yoda_histogram2, ...}}}
        Usually only "0" exists as the innermost key.
    refhsitos : dict
        Dictionary of the reference analysis data YODA histograms.
    filelist : list[str]

    filenames : list[str]
    
    plotoptions : dict[str, dict[str, str]]
    
    Returns
    -------
    outputdict : dict
        
    """
    outputdict = {'rivet': yamlparser.get_plot_configs(plot_id, plotdirs=plotdirs, config_files=config_files)}
    # TODO: Will there ever be preexisting histograms?
    outputdict['histograms'] = {}
    if plot_id in refhistos:
        with io.StringIO() as filelike_str:
            yoda.writeFLAT(refhistos[plot_id], filelike_str)
            outputdict['histograms'][reftitle] = {histogram_str_name: literal(filelike_str.getvalue())}
            # TODO: add plot options here as well.
            # TODO: refactor so that same function is applied to refhistos and mchistos

    for filename, mchistos_in_file in mchistos.items():
        outputdict['histograms'][filename] = {}
        # TODO: will there ever be multiple histograms with same ID here? Rewrite _get_histos?
        for histogram in mchistos_in_file[plot_id].values():
            # TODO: Probably exists a more efficient way of doing this. Just looping over all settings maybe.
            outputdict['histograms'][filename].update(plotoptions.get(filename, {}))
            outputdict['histograms'][filename].update(plotoptions.get('PLOT', {}))

            with io.StringIO() as filelike_str:
                yoda.writeFLAT(histogram, filelike_str)
                # TODO: Check with rivet-cmphistos that the name change is correct
                outputdict['histograms'][filename][histogram_str_name] = literal(filelike_str.getvalue())
    return outputdict


def mkoutdir(outdir):
    """Function to make output directories"""
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except:
            msg = "Can't make output directory '%s'" % outdir
            raise Exception(msg)
    if not os.access(outdir, os.W_OK):
        msg = "Can't write to output directory '%s'" % outdir
        raise Exception(msg)


def _write_output(output, h, hier_output, outdir):
    "Choose output file name and dir"
    if hier_output:
        hparts = h.strip("/").split("/", 1)
        ana = "_".join(hparts[:-1]) if len(hparts) > 1 else "ANALYSIS"
        outdir = os.path.join(outdir, ana)
        outfile = '%s.yaml' % hparts[-1].replace("/", "_")
    else:
        hparts = h.strip("/").split("/")
        outfile = '%s.yaml' % "_".join(hparts)
    mkoutdir(outdir)
    outfilepath = os.path.join(outdir, outfile)
    with open(outfilepath, 'w') as yaml_file:
        yaml.dump(output, yaml_file, indent=4,  default_flow_style=False)

        
def make_yamlfiles(args, path_pwd=True, reftitle='Data', 
                   rivetrefs=True, path_patterns=(), 
                   path_unpatterns=(), plotinfodirs=[], 
                   style='default', plot_features='',
                   config_files=[], hier_output=False, outdir='.',
                   rivetplotpaths=True, rc_params={}
                  ):
    """
    Create .yaml files that can be parsed by rivet-make-plot
    Each output .yaml file corresponds to one analysis which contains all MC histograms and a reference data histogram.
    Warning: still in development.
    
    Parameters
    ----------
    args : Iterable[str]
        Non-keyword arguments that were previously passed to rivet-cmphistos. 
        E.g., ['mc1.yoda', 'mc2.yoda:Title=example title'] 
        TODO: change this input to filelist, filenames, plotoptions instead?
    path_pwd : bool
        Search for plot files and reference data files in current directory.
    reftitle : str
        Legend name of the reference data in the plots.
    rivetrefs : bool
        If False, don't use Rivet reference data files
    path_patterns : Iterable[str]
        Only write out histograms whose $path/$name string matches these regexes.
        The argument may also be a text file.
    path_unpatterns : Iterable[str]
        Exclude histograms whose $path/$name string matches these regexes
    TODO: path_patterns, path_unpatterns have probably not been implemented yet.
    plotinfodirs : list[str]
        Directory which may contain plot header information (in addition to standard Rivet search paths).
    style : str TODO
        Set the style of all plots. 
        Either a .yaml file name or a name of a builtin style (e.g. 'default')
    plot_features : str TODO
        Settings that will be included in the "plot features" section of the YAML file.
        The string has the format "key=value:key2=value2"..., i.e., a similar format as the strings in args.
    config_files : list[str]
        Additional plot config file(s). 
        Settings will be included in the output configuration. 
        ~/.make-plots will automatically be added.
    hier_output : bool
        Write output .yaml files into a directory hierarchy which matches the analysis paths.
    outdir : str
        Write yaml files into this directory.
    rivetplotpaths : bool
        Search for .plot files in the standard Rivet plot paths.
    rc_params : dict[str, str] TODO
        Additional rc params added to all output .yaml files. 
    Returns
    -------
    None
    
    Raises
    ------
    IOError
        If the program does not have read access to .plot or .yoda files, or if it cannot write the output .yaml files.
    """
    # Code from rivet-cmphistos >>> 
    # TODO: clean and refactor rivet-cmphistos code

    ## Add pwd to search paths
    if path_pwd:
        rivet.addAnalysisLibPath(os.path.abspath("."))
        rivet.addAnalysisDataPath(os.path.abspath("."))
    
    # Add .make-plots which contains extra plotting configurations.
    config_files.append('~/.make-plots')
    ## Split the input file names and the associated plotting options
    ## given on the command line into two separate lists
    filelist, filenames, plotoptions = _parse_args(args)
    
    ## Check that the files exist
    for f in filelist:
        if not os.access(f, os.R_OK):
            raise IOError("Error: cannot read from %s" % f)
    
    plotdirs = plotinfodirs + [os.path.abspath(os.path.dirname(f)) for f in filelist] + (rivet.getAnalysisPlotPaths() if rivetplotpaths else [])

    # Create a list of all histograms to be plotted, and identify if they are 2D histos (which need special plotting)
    # TODO: implement 2D histogram special settings. 
    refhistos, mchistos = _get_histos(filelist, filenames, plotoptions, path_patterns, path_unpatterns)
    
    hpaths = []
    for aos in mchistos.values():
        for p in aos.keys():
            ps = rivet.stripOptions(p)
            if ps and ps not in hpaths:
                hpaths.append(ps)
            # Only use the first histogram
            firstaop = aos[p][sorted(aos[p].keys())[0]]

    # Unique list of analyses
    anas = list(set([x.split("/")[1] for x in hpaths]))

    ## Take reference data from the Rivet search paths, if there is not already
    if rivetrefs:
        refhistos2 = _get_rivet_ref_data(anas, path_patterns, path_unpatterns)
        refhistos2.update(refhistos)
        refhistos = refhistos2

    ## Purge unmatched ref data entries to save memory
    keylist = list(refhistos.keys())
    for refhpath in keylist:
        if refhpath not in hpaths:
            del refhistos[refhpath]
    # <<< end of code from rivet-cmphistos
    
    # Write each file 
    for plot_id in hpaths:
        outputdict = _make_output(plot_id, plotdirs, config_files, mchistos, refhistos, reftitle, filelist, filenames, plotoptions)
        
        ## Make the output and write to file
        _write_output(outputdict, plot_id, hier_output=hier_output, outdir=outdir)
