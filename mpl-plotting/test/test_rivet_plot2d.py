import numpy as np
import pytest
import yoda

from rivet_plot import rivet_plot
from rivet_plot2d import plot_2Dhist

# This is used instead of relative paths so that pytest can be run from anywhere.
import os
io_file_dir = os.path.join(os.path.dirname(__file__), '2d-io-files/')

# TODO more tests for different settings.
def default_hist_features(size):
    return tuple([{'Title': 'Data'}]+[{'Title': 'mc' + str(i)} for i in range(1, size)])

def empty_hists(size):
    return tuple(yoda.Histo2D(1, -1, 1, 1, -1, 1).mkScatter() for i in range(size))

def filled_hists(size, nevents=1000, histo_args=(30, -1, 2, 20, -1, 1)):
    hists = []
    np.random.seed(1)

    for data in np.stack((np.random.normal(np.mean(histo_args[1:2]), 1, size=(size, nevents)), 
                          np.random.normal(np.mean(histo_args[4:5]), 1, size=(size, nevents))), 
                        axis=-1):
        h = yoda.Histo2D(*histo_args)
        for event in data:
            h.fill(event[0], event[1])
        hists.append(h.mkScatter())
    return tuple(hists)

def empty_dicts(size):
    return tuple({} for i in range(size))


# 3-empty-hists gives strange colorbar. This is ok though, since people should not plot empty histograms anyway
@pytest.mark.parametrize(
    'hist_data,hist_features,plot_features,filename',
    [
        (filled_hists(3), default_hist_features(3), {'RatioPlot': True, 'LogZ': False}, '3-filled-hists'),
        (filled_hists(3), default_hist_features(3), {'ShowZero': False, 'LogZ': False}, '3-filled-hists-no-showzero'),
        (filled_hists(4), default_hist_features(4), {'LogZ': False}, '4-filled-hists'),
        (
            filled_hists(3, 100, histo_args=(5, 1, 2, 6, 0.09, 1.1)), default_hist_features(3), 
            {'2DType': 'heatmap', '2DIndividual': False, 'LogX': True, 'LogY': True, 'LogZ': True}, 
            '3-filled-hists-log'
        ),
        (empty_hists(3), default_hist_features(3), {'LogZ': False}, '3-empty-hists'),
    ]
)
def test_no_error_plot_2Dhist(hist_data, hist_features, plot_features, filename):
    yaml_dict = {'plot features': plot_features}
    plot_2Dhist(hist_data, hist_features, yaml_dict, io_file_dir + filename)


@pytest.mark.parametrize(
    'hist_data,hist_features,plot_features,error',
    [
        (empty_hists(3), default_hist_features(3), {'2DType': 'heatmap-something'}, ValueError),
        (empty_hists(3), default_hist_features(3), {'2DType': ''}, ValueError),
        (empty_hists(3), default_hist_features(3), {'2DType': 'contour'}, ValueError),
        (empty_hists(3), empty_dicts(3), {'2DType': 'heatmap', 'LogZ': False}, KeyError)
    ]
)
def test_error_rivet_plot2d(hist_data, hist_features, plot_features, error):
    yaml_dict = {'plot features': plot_features}
    with pytest.raises(error):
        plot_2Dhist(hist_data, hist_features, yaml_dict, 'name should not be important')

# TODO test other file formats than png

@pytest.mark.parametrize(
    'filename',
    [
        io_file_dir + 'log-axis.yaml', 
        io_file_dir + 'only-ref.yaml',
        io_file_dir + 'all-default.yaml',
        io_file_dir + 'surface.yaml',
        io_file_dir + 'one-figure.yaml',
        io_file_dir + 'one-figure-custom-cmap.yaml',
        io_file_dir + 'one-figure-surface.yaml',
        io_file_dir + 'one-figure-surface-view.yaml'
    ]
)
def test_no_error_rivet_plot(filename):
    """Test that the rivet_plot code runs without errors for 2D histograms.

    Parameters
    ----------
    filename : str
        File name of a yaml file that will be plotted.
    
    Raises
    ------
    Exception
        If something goes wrong. Images need to be reviewed manually.
    
    Notes
    -----
    This test function (or a variant) can later be used in conjuction with pytest-mpl
    """
    # outputdir is '/' because file names are absolute paths.
    rivet_plot(filename, filename.rstrip('.yaml'), outputdir='/')


# To profile code, run: 
#  python -m cProfile -o rivet_plot.pstat test_rivet_plot2d.py