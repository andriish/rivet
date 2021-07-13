from __future__ import print_function
import os, re
from ruamel import yaml
import rivet
import constants
from old_plotfile_converter import parse_old_plotfile


# This class, function and call to add_representer makes it so that all objects with type Literal will be printed to a .yaml file as a string block. 
class Literal(str):
    """A small wrapper class used to print histograms as multiline strings in a .yaml file."""
    pass

def literal_presenter(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(Literal, literal_presenter)


def _parse_yaml_plotfile(filename, hpath):
    """Parse a .plot file with the yaml format and return the settings.

    Parameters
    ----------
    filename : str
        The name of a .plot file with yaml-like syntax.
    hpath : str
        The histogram path, usually with the format /AnalysisID/HistogramID .
    
    Returns
    -------
    plotfile_configs : dict
        All the settings from the file
    """
    plotfile_configs = {}
    with open(filename, 'r') as file:
        for configs in yaml.safe_load_all(file):
            if re.match(configs[constants.name_key], hpath):   # TODO: old code uses match instead of fullmatch. Is this intentional?
                plotfile_configs.update(configs[constants.plot_setting_key])
    return plotfile_configs


def _is_old_format(filename):
    """Check if a file contains the old .plot format or the new yaml format.

    Parameters
    ----------
    filename : str
        Path to the .plot file that will be checked.
    
    Returns
    -------
    is_old : bool
        True if the file is the old format. False otherwise.

    Notes
    -----
    This does not make an explicit check that the format is correct yaml syntax.
    Instead, it only checks whether the file contains #BEGIN PLOT or not.
    TODO: optimize if needed. Current implementation is probably slow.
    """
    is_old = False
    with open(filename) as file:
        for line in file:
            if re.match(r'^(#*\s*)?BEGIN (\w+) ?(\S+)? ?(\w+)?', line):
                is_old = True
                break
    return is_old


def _get_matching_plot_configs_from_file(hpath, plotfilepath):
    """Open plotfilepath and get the settings with the corresponding hpath ID.
    
    Parameters
    ----------
    hpath : str
        The histogram path, usually with the format /AnalysisID/HistogramID.
    plotfilepath : str
        The path to a .plot file, either with yaml syntax or in the old format.

    Returns
    -------
    plot_configs : dict
        Dictionary of settings for the corresponding histogram from the .plot file.
        Will return an empty dict if the hpath could not be found in the file.
    """
    plot_configs = {}
    if not os.access(plotfilepath, os.R_OK):
        return {}
        
    if _is_old_format(plotfilepath):
        new_plot_settings = parse_old_plotfile(hpath, plotfilepath)  # BUG: Were these accidentally switched?
    else:
        new_plot_settings = _parse_yaml_plotfile(plotfilepath, hpath)
    plot_configs.update(new_plot_settings)

    return plot_configs


def get_plot_configs(hpath, plotdirs=[], config_files=[]):
    """Get all settings for the hpath analysis by reading through the settings of plotdirs and config_files
    
    Parameters
    ----------
    hpath : str
        The histogram path,  usually with the format /AnalysisID/HistogramID .
    plotdirs : list[str]
        Directories containing .plot files. The settings in the files with name AnalysisID.plot will be parsed and added to plot_configs.
    config_files : Iterable[str]
        Extra .plot files with settings. If there are sections in these yaml files with name /AnalysisID/HistogramID, these settings will be added to plot_configs.
    
    Returns
    -------
    plot_configs : dict
        Dict containing all the settings from the .plot files that match the hpath ID. 
    
    Note
    ----
    Warning: if a key exists in multiple configs, the last setting will included in plot_configs. 
    As an example, is multiple .plot files have the `xlabel` setting, the xlabel specified in the last file will be used. 
    I hope this is an expected behavior. TODO: ask mentors
    """
    # Remove duplicates
    plotdirs = list(set(plotdirs))
    
    id_parts = rivet.AOPath(hpath).basepathparts()

    plot_configs = {}
    for plotdir in plotdirs:
        plotfilepath = os.path.join(plotdir, id_parts[0] + constants.file_extension)
        plotfile_configs = _get_matching_plot_configs_from_file(hpath, plotfilepath)   # TODO: Can I pass hpath here or will that cause problems with /REF?
        plot_configs.update(plotfile_configs)
        
    for plotfilepath in config_files:
        plotfile_configs = _get_matching_plot_configs_from_file(hpath, plotfilepath)   # TODO: Can I pass hpath here or will that cause problems with /REF?
        plot_configs.update(plotfile_configs)
        
    return plot_configs


def _mkoutdir(outdir):
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


def write_output(output, h, hier_output, outdir):
    "Choose output file name and dir"
    if hier_output:
        hparts = h.strip("/").split("/", 1)
        ana = "_".join(hparts[:-1]) if len(hparts) > 1 else "ANALYSIS"
        outdir = os.path.join(outdir, ana)
        outfile = '%s.yaml' % hparts[-1].replace("/", "_")
    else:
        hparts = h.strip("/").split("/")
        outfile = '%s.yaml' % "_".join(hparts)
    _mkoutdir(outdir)
    outfilepath = os.path.join(outdir, outfile)
    with open(outfilepath, 'w') as yaml_file:
        yaml.dump(output, yaml_file, indent=4,  default_flow_style=False)


def read_yamlfile(filename):
    """Read a .yaml file with a single section inside it (i.e., no multiple sections).

    Parameters
    ----------
    filename : 
        Name of yaml file.
    
    Returns
    -------
    content : Any
        Content of the yaml file.

    Note
    ----
    This is a thin wrapper around the yaml backend. 
    This function exists so that the yaml-reading backend can be switched (e.g., if one backend does not support python 2).
    It might be removed in the future.
    """
    with open(filename) as file:
        content = yaml.safe_load(file)
    return content

