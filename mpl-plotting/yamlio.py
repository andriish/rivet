from __future__ import print_function
import os, re
from ruamel import yaml
import rivet

# TODO: remove these variables once the names have been properly decided
name_key = 'name'
plot_setting_key = 'plot features'
file_extension = '.plot'

# This class, function and call to add_representer makes it so that all objects with type literal will be printed to a .yaml file as a string block. 
class literal(str):
    """A small wrapper class used to print histograms as multiline strings in a .yaml file."""
    pass

def literal_presenter(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(literal, literal_presenter)


def _get_matching_plot_configs_from_file(hpath, plotfilepath): # TODO: better variable names
    """Open plotfilepath, which is a .plot file with yaml syntax, parse it and get the settings with the hpath ID.
    """
    plot_configs = {}
    if not os.access(plotfilepath, os.R_OK):
        return {}
        
    # TODO: do some preprocessing (or processing at the same time) with variables.
    with open(plotfilepath, 'r') as plot_config_file:
        for configs in yaml.safe_load_all(plot_config_file):
            if hpath == configs[name_key]:   # TODO: change to regex later
                plot_configs.update(configs[plot_setting_key])
    return plot_configs


# TODO: turn this into a class?
def get_plot_configs(hpath, plotdirs=[], config_files=[]):
    """Get all settings for the hpath analysis by reading through the settings of plotdirs and config_files
    
    Parameters
    ----------
    hpath : str
        The histogram path, with format /AnalysisID/HistogramID .
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
        plotfilepath = os.path.join(plotdir, id_parts[0] + file_extension)
        plotfile_configs = _get_matching_plot_configs_from_file(hpath, plotfilepath)   # TODO: Can I pass hpath here or will that cause problems with /REF?
        plot_configs.update(plotfile_configs)
        
    for plotfilepath in config_files:
        plotfile_configs = _get_matching_plot_configs_from_file(hpath, plotfilepath)   # TODO: Can I pass hpath here or will that cause problems with /REF?
        plot_configs.update(plotfile_configs)
        
    return plot_configs


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
    mkoutdir(outdir)
    outfilepath = os.path.join(outdir, outfile)
    with open(outfilepath, 'w') as yaml_file:
        yaml.dump(output, yaml_file, indent=4,  default_flow_style=False)


# Test code
if __name__ == '__main__':
    # TODO: plotdirs will be changed to rivet.getAnalysisPlotPaths() once .plot files have been put there
    print(get_plot_configs('/ALICE_2010_S8625980/d03-x01-y01', [os.getcwd()]))
