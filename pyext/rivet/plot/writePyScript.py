# TODO

import os
import yoda

from rivet.plot.mpl_tools.yamlio import read_yamlfile

def testPrint():
  print("hi")

def writePyScript1D():
  # TODO
  pass

def writePyScript2D():
  #TODO 
  pass

def writePyScript(yaml_dict, plot_name, outdir):
  """ 
  write a python script that only relies on YODA
  functionalities
  """
    
  outPyName = os.path.join(outdir, plot_name.strip('/')) + '.py'

  # contents to write out to .py script
  mplScript = """#! /usr/bin/env python"""

  mplScript +='''
"""
This script relies on Python and allows you to recreate
the plot made with rivet-mkhtml. However now with the flexibility
to change the plot features, for example log/lin axes,
line colours or the position of the legend.
"""
'''

  # add imports
  mplScript += """
import matplotlib.pyplot as plt
from yoda import mkCanvas,mkPlots
"""
  
  # TODO 0
  # distinguish between making a 1D or 2D plot
  # see how rivet_plot() currently does this
  # for example:
  # if <1D>: mplScript   += writePy1D(...)
  # elif <2D>: mplScript += writePy2D(...)
   
  # TODO 1 (inside writePy1D function)
  # get numerical values of .yoda files and reference data
  # save this in some kind of lists. Can't be in-memory YODA objects
  # but need to be converted to ascii to live in the .py script
  # save as dictionary, e.g.
  # plotData = {'refData' : <list>, 'yoda1' : <list>, 'yoda2' : <list>}
  # or use pickle here?
  
  with open(outPyName, "a") as f:
    f.write(mplScript)
    #f.write("import matplotlib.pyplot as plt")

  return(outPyName)
   
# ultimately, have something produced of this form:
  """
    #!/usr/bin/python
    import matplotlib.pyplot as plt
    from YODA import mkplots, mkCanvas

    style = 'rivet'
    fig,ax = yoda.mkCanvas(style, ratio='True') # --> this makes the canvas, tick size, etc. etc.
    
    ax.set_ylim(5, 20)
    ax.set_yscale('log')
    ax.set_xscale('lin')
    plt.rcParams['legend.loc'] = "upper right"

    # these lists need to contain explicit numerical values, otherwise can't
    # run it separately from rivet
    plotData = {'refData' : <list>, 'yoda1' : <list>, 'yoda2' : <list>}
    bins = <list>
    yoda.mkPlots(fig, ax, bins, plotData['refData'], color='black', type='Data')
    yoda.mkPlots(fig, ax, bins, plotData['yoda1'],    color='green', type='MC') 
    yoda.mkPlots(fig, ax, bins, plotData['yoda2'],    color='green', type='MC') 
    
    VECTORFORMAT = 'pdf' # command line option with rivet-mkhtml
    plt.savefig(f'/path/to/file.{VECTORFORMAT}')
  """

  
  pass  
