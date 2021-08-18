"""Module creates a rivet-style plot as a pdf."""
import os
import io

import yoda
from yamlio import read_yamlfile
from rivet_plot1d import plot_1Dhist
from rivet_plot2d import plot_2Dhist


def rivet_plot(yaml_file, plot_name, outputdir='.'):
    """Create plot from the yaml file in the Rivet style.

    Parameters
    ----------
    yaml_file : str or dict
        Either the file name of the yaml file or a dictionary object containing the plot data.
    plot_name : str
        Name of the created plot file.
    outputdir : str, optional
        Name of the relative output directory. Default is the current directory.
    """
    # Parse yaml file for histogram data
    if isinstance(yaml_file, str):  # If the file name of the yaml file is passed
        yaml_dicts = read_yamlfile(yaml_file)
        hist_data = _parse_yoda_hist(yaml_dicts)
    else:  # If the dictionary object is passed
        yaml_dicts = yaml_file
        hist_data = [
            hist_dict['yoda'].mkScatter()
            if not isinstance(hist_dict, (yoda.Scatter1D, yoda.Scatter2D, yoda.Scatter3D)) else hist_dict['yoda']
            for hist_dict in yaml_dicts['histograms'].values()
        ]

    hist_features = [val for val in yaml_dicts['histograms'].values()]
    output_filename = os.path.join(outputdir, plot_name.strip('/'))

    # Ensure reference histogram is first in list since dicts are not ordered in Python 2, not needed for Python 3
    for i, h in enumerate(hist_features):
        if h.get('IsRef'):
            if i != 0:  # Not necessary since the reference histogram is already the first element
                hist_data.insert(0, hist_data.pop(i))
                hist_features.insert(0, hist_features.pop(i))
            break

    if all(isinstance(h, yoda.Scatter1D) for h in hist_data) or all(isinstance(h, yoda.Scatter2D) for h in hist_data):
        plot_1Dhist(hist_data, hist_features, yaml_dicts, output_filename)
    elif all(isinstance(h, yoda.Scatter3D) for h in hist_data):
        plot_2Dhist(hist_data, hist_features, yaml_dicts, output_filename)
    else:
        print('Error with Class types:', [type(h) for h in hist_data])
        raise NotImplementedError('Class type cannot be plotted yet')


def _parse_yoda_hist(yaml_dicts):
    """Read yoda string and return yoda object."""
    hist_data = []
    for hist_dict in yaml_dicts['histograms'].values():
        with io.StringIO(hist_dict['yoda']) as file_like:
            hist_data.append(yoda.readYODA(
                file_like, asdict=False)[0].mkScatter())
    return hist_data
