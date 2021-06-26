import rivet, yoda
import os, glob, io
import yamlparser
from __future__ import print_function

# This makes it so that all objects with type literal will be printed to a .yaml file as a string block. 

class literal(str):
    """A small wrapper class used to print histograms as string blocks in a .yaml file."""
    pass

def literal_presenter(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(literal, literal_presenter)


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
def _get_histos(filelist, filenames):
    """Loop over all input files. Only use the first occurrence of any REF-histogram
    and the first occurrence in each MC file for every MC-histogram."""
    refhistos, mchistos = {}, {}
    for infile, inname in zip(filelist, filenames):
        mchistos.setdefault(inname, {})
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
            defaultWeightName = getDefaultVariation(plotoptions[inname])
            if aop.isref() and basepath not in refhistos:
                ao.setPath(aop.varpath(keepref=False, defaultvarid=defaultWeightName))
                refhistos[basepath] = ao
            else: #if basepath not in mchistos[infile]:
                mchistos[inname].setdefault(basepath, {})[aop.varid(defaultWeightName)] = ao

    return refhistos, mchistos


# GSoC TODO: keep
def _get_rivet_ref_data(anas=None):
    "Find all Rivet reference data files"
    refhistos = {}
    rivet_data_dirs = rivet.getAnalysisRefPaths()
    dirlist = []
    for d in rivet_data_dirs:
        import glob
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
            analysisobjects = yoda.read(infile, patterns=args.PATHPATTERNS, unpatterns=args.PATHUNPATTERNS)
            for path, ao in analysisobjects.items():
                aop = rivet.AOPath(ao.path())
                if aop.isref():
                    ao.setPath(aop.basepath(keepref=False))
                    refhistos[ao.path()] = ao
    return refhistos


# GSoC TODO: keep
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


# GSoC TODO: change to new dat, py or yaml format. Currently yaml
def _write_output(output, h):
    "Choose output file name and dir"
    if args.HIER_OUTPUT:
        hparts = h.strip("/").split("/", 1)
        ana = "_".join(hparts[:-1]) if len(hparts) > 1 else "ANALYSIS"
        outdir = os.path.join(args.OUTDIR, ana)
        outfile = '%s.dat' % hparts[-1].replace("/", "_")
    else:
        hparts = h.strip("/").split("/")
        outdir = args.OUTDIR
        outfile = '%s.yaml' % "_".join(hparts)
    mkoutdir(outdir)
    outfilepath = os.path.join(outdir, outfile)
    with open(outfilepath, 'w') as yaml_file:
        yaml.dump(output, yaml_file, indent=4)

        
def make_yamlfiles(mcfilenames, *args, **kwargs):
    """
    Create .yaml files that can be parsed by rivet-make-plot
    Each output .yaml file corresponds to one analysis which contains all MC histograms and a reference data histogram.
    Warning: still in development.
    
    Parameters
    ----------
    mcfilenames : list[str]
        yoda files with MC histograms inside them. 
    args, kwargs : 
        options that will be added in the future, such as style (e.g., ATLAS), extra plotting options, hierout option, and extra files with plotting options.
        
    Returns
    -------
    None
    """
    # Code from rivet-cmphistos >>> 
    # TODO: clean and refactor rivet-cmphistos code
    
    ## Add pwd to search paths
    if args.PATH_PWD:
        rivet.addAnalysisLibPath(os.path.abspath("."))
        rivet.addAnalysisDataPath(os.path.abspath("."))

    ## Split the input file names and the associated plotting options
    ## given on the command line into two separate lists
    # GSoC TODO: keep as is for now
    filelist, filenames, plotoptions = _parse_args(mcfilenames)

    ## Check that the files exist
    for f in filelist:
        if not os.access(f, os.R_OK):
            print("Error: cannot read from %s" % f)
            sys.exit(1)
    
    # TODO: add rivet.getAnalysisPath
    plotdirs = args.PLOTINFODIRS + [os.path.abspath(os.path.dirname(f)) for f in filelist] + rivet.getAnalysisPlotPaths()

    ## Create a list of all histograms to be plotted, and identify if they are 2D histos (which need special plotting)
    try:
        refhistos, mchistos = _get_histos(filelist, filenames)
    except IOError as e:
        print("File reading error: ", e.strerror)
        exit(1)
    hpaths, h2ds = [], []
    for aos in mchistos.values():
        for p in aos.keys():
            ps = rivet.stripOptions(p)
            if ps and ps not in hpaths:
                hpaths.append(ps)
            firstaop = aos[p][sorted(aos[p].keys())[0]]
            # TODO: Would be nicer to test via isHisto and dim or similar, or yoda.Scatter/Histo/Profile base classes
            if type(firstaop) in (yoda.Histo2D, yoda.Profile2D, yoda.Scatter3D) and ps not in h2ds:
                h2ds.append(ps)

    # Unique list of analyses
    anas = list(set([x.split("/")[1] for x in hpaths]))
    #print (anas)

    ## Take reference data from the Rivet search paths, if there is not already
    if args.RIVETREFS:
        try:
            refhistos2 = _get_rivet_ref_data(anas)
        except IOError as e:
            print("File reading error: ", e.strerror)
            exit(1)
        refhistos2.update(refhistos)
        refhistos = refhistos2

    ## Purge unmatched ref data entries to save memory
    for refhpath in refhistos:
        if refhpath not in hpaths:
            del refhistos[refhpath]
    # <<< end of code from rivet-cmphistos
    
    # Write each file 
    for analysis_id in hpaths:
        # TODO: move to separate function?
        outputdict = yamlparser.get_plot_configs(analysis_id, plotdirs=plotdirs, extra_files=extra_files) # TODO: add extra_files and plotdirs from args
        # Will there ever be preexisting histograms?
        outputdict['histograms'] = {}
        # TODO: add option to rename Data here
        outputdict['histograms']['Data'] = refhistos[analysis_id]
        
        for filename, mchistos_in_file in mchistos.items():
            # TODO: move this section to separate function
            # TODO: will there ever be multiple histograms with same ID here? Rewrite _get_histos?
            for histogram in mchistos_in_file[analysis_id].values():
                sio = io.StringIO()
                yoda.writeFLAT(histogram, sio)
                # TODO: add option to change name of histogram here
                outputdict['histograms'][filename.strip('.yoda')] = sio.getvalue()
        
        ## Make the output and write to file
        _write_output(outputdict, hpath)

    