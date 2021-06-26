from __future__ import print_function
import rivet, yoda
import os, glob, io
import yamlparser 
# TODO: move all yaml-related things to a yamlio.py file so that if parser needs to be changed, it is only in that file?
from ruamel import yaml

# TODO: replace exits with raise exception instead since this is a function, not a command line tool.

# This makes it so that all objects with type literal will be printed to a .yaml file as a string block. 

class literal(str):
    """A small wrapper class used to print histograms as string blocks in a .yaml file."""
    pass

def literal_presenter(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(literal, literal_presenter)


# GSoC TODO: keep
def sanitiseString(s):
    #s = s.replace('_','\\_')
    #s = s.replace('^','\\^{}')
    #s = s.replace('$','\\$')
    s = s.replace('#','\\#')
    s = s.replace('%','\\%')
    return s


# GSoC TODO: keep
def getDefaultVariation(options):
    name = "0"
    for val in getFileOptions(options, "DefaultWeight").values():
        name = val
    return name


# GSoC TODO: keep
def getFileOptions(options, tags):
    ret = { }
    if 'list' not in str(type(tags)):    # GSoC TODO: Why is this not just `if isinstance(tags, list):`?
        tags = [ tags ]
    for opt in options:
        key, val = opt.split('=', 1)
        if key in tags:
            ret[key] = val
    return ret


# GSoC TODO: one should be able to specify plot options here but as matplotlib args.
#     Old plot options should be converted to matplotlib plot options.
#     These might just be kwargs passed to the underlying plotting function.  
#     As before, the `PLOT` argument specifies args for all yoda files.
#     Example: `rivet-cmphistos mc1.yoda:label='first monte carlo generator' mc2.yoda PLOT:linestyle='--'`
# Problem: some matplotlib line styles contain ':', which would not work with current code. Might have to rewrite.
def _parse_args(args):
    """Look at the argument list and split it at colons, in order to separate
    the file names from the plotting options. Store the file names and
    file specific plotting options."""
    filelist = []
    filenames = []
    plotoptions = {}
    for a in args:
        asplit = a.split(':')
        path = asplit[0]
        if path != "PLOT":
            filelist.append(path)
            filenames.append(path)
        plotoptions[path] = []
        has_title = False
        has_name = ""
        for i in range(1, len(asplit)):
            ## Add 'Title' if there is no = sign before math mode
            if '=' not in asplit[i] or ('$' in asplit[i] and asplit[i].index('$') < asplit[i].index('=')):
                asplit[i] = 'Title=%s' % asplit[i]
            if asplit[i].startswith('Title='):
                has_title = True
            plotoptions[path].append(asplit[i])
            if asplit[i].startswith('Name=') and path != "PLOT":
                has_name = asplit[i].split('=', 1)[1]
                filenames[-1] = has_name
        if has_name != "":
            plotoptions[has_name] = plotoptions[path]
            del plotoptions[path]
        if path != "PLOT" and not has_title:
            plotoptions[has_name if has_name != "" else path].append('Title=%s' % sanitiseString(os.path.basename( os.path.splitext(path)[0] )) )
    return filelist, filenames, plotoptions


# GSoC TODO: keep as is for now. Improve with better "regex" support if time permits. Maybe this kind of regex should only be a "feature" for the API? 
def _get_histos(filelist, filenames, plotoptions={}, path_patterns=(), path_unpatterns=()):
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
            defaultWeightName = getDefaultVariation(plotoptions[inname])
            if aop.isref() and basepath not in refhistos:
                ao.setPath(aop.varpath(keepref=False, defaultvarid=defaultWeightName))
                refhistos[basepath] = ao
            else: #if basepath not in mchistos[infile]:
                mchistos[inname].setdefault(basepath, {})[aop.varid(defaultWeightName)] = ao

    return refhistos, mchistos


def _get_rivet_ref_data(anas=None, path_patterns=(), path_unpatterns=()):
    "Find all Rivet reference data files"
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


# TODO docstr 
def _write_output(output, h, hier_output=False, outdir='rivet-plots'):
    "Choose output file name and dir"
    if hier_output:
        hparts = h.strip("/").split("/", 1)
        ana = "_".join(hparts[:-1]) if len(hparts) > 1 else "ANALYSIS"
        outdir = os.path.join(outdir, ana)
        outfile = '%s.dat' % hparts[-1].replace("/", "_")
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
                   config_files=[], hier_output=False,
                   rivetplotpaths=True,
                   **kwargs):
    """
    Create .yaml files that can be parsed by rivet-make-plot
    Each output .yaml file corresponds to one analysis which contains all MC histograms and a reference data histogram.
    Warning: still in development.
    
    Parameters
    ----------
    args : list[str]
        Non-keyword arguments that were previously passed to plot. 
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
    plotinfodirs : list[str]
        Directory which may contain plot header information (in addition to standard Rivet search paths).
    config_files : list[str]
        Additional plot config file(s). Settings will be included in the output configuration. ~/.make-plots will automatically be added.
    hier_output : bool
        Write output dat files into a directory hierarchy which matches the analysis paths.
    rivetplotpaths : bool
    kwargs : 
        options that will be added in the future, such as style (e.g., ATLAS), extra plotting options, hierout option, and extra files with plotting options.
        
    Returns
    -------
    None
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
    # TODO: keep as is for now. Add support for plotoptions
    filelist, filenames, plotoptions = _parse_args(args)
    
    ## Check that the files exist
    for f in filelist:
        if not os.access(f, os.R_OK):
            raise IOError("Error: cannot read from %s" % f)
    
    # TODO: add rivet.getAnalysisPath?
    plotdirs = plotinfodirs + [os.path.abspath(os.path.dirname(f)) for f in filelist] + (rivet.getAnalysisPlotPaths() if rivetplotpaths else [])

    ## Create a list of all histograms to be plotted, and identify if they are 2D histos (which need special plotting)
    try:
        refhistos, mchistos = _get_histos(filelist, filenames, plotoptions)
    except IOError as e:
        print("File reading error: ", e.strerror)
        exit(1)
        
    # h2ds is currently not used
    hpaths, h2ds = [], []
    for aos in mchistos.values():
        for p in aos.keys():
            ps = rivet.stripOptions(p)
            if ps and ps not in hpaths:
                hpaths.append(ps)
            # Only use the first histogram
            firstaop = aos[p][sorted(aos[p].keys())[0]]
            # TODO: Would be nicer to test via isHisto and dim or similar, or yoda.Scatter/Histo/Profile base classes
            if type(firstaop) in (yoda.Histo2D, yoda.Profile2D, yoda.Scatter3D) and ps not in h2ds:
                h2ds.append(ps)

    # Unique list of analyses
    anas = list(set([x.split("/")[1] for x in hpaths]))

    ## Take reference data from the Rivet search paths, if there is not already
    if rivetrefs:
        try:
            refhistos2 = _get_rivet_ref_data(anas)
        except IOError as e:
            print("File reading error: ", e.strerror)
            exit(1)
        refhistos2.update(refhistos)
        refhistos = refhistos2

    ## Purge unmatched ref data entries to save memory
    keylist = list(refhistos.keys())
    for refhpath in keylist:
        if refhpath not in hpaths:
            del refhistos[refhpath]
    # <<< end of code from rivet-cmphistos
    
    # Write each file 
    # TODO: rename plot_id
    for plot_id in hpaths:
        # TODO: move to separate function?
        outputdict = {'rivet': yamlparser.get_plot_configs(plot_id, plotdirs=plotdirs, config_files=config_files)} # TODO: add config_files and plotdirs from args
        # Will there ever be preexisting histograms?
        outputdict['histograms'] = {}
        # TODO: add option to rename Data here
        if plot_id in refhistos:
            with io.StringIO() as filelike_str:
                yoda.writeFLAT(refhistos[plot_id], filelike_str)
                outputdict['histograms'][reftitle] = literal(filelike_str.getvalue())
        
        for filename, mchistos_in_file in mchistos.items():
            # TODO: move this section to separate function
            # TODO: will there ever be multiple histograms with same ID here? Rewrite _get_histos?
            for histogram in mchistos_in_file[plot_id].values():
                with io.StringIO() as filelike_str:
                    yoda.writeFLAT(histogram, filelike_str)
                    # TODO: change name of histogram here
                    outputdict['histograms'][filename] = literal(filelike_str.getvalue())

        ## Make the output and write to file
        _write_output(outputdict, plot_id)


# Test code
if __name__ == '__main__':
    make_yamlfiles(['mc1.yoda', 'mc2.yoda'], hier_out=True, rivetplotpaths=False)