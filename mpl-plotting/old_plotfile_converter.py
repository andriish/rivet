from __future__ import print_function
import os, re, logging
from rivet.aopaths import AOPath
from rivet.util import texpand
# TODO: make it a class, as it was originally, instead?

pat_begin_block = re.compile(r'^(#*\s*)?BEGIN (\w+) ?(\S+)?')
pat_begin_name_block = re.compile(r'^(#*\s*)?BEGIN (\w+) ?(\S+)? ?(\w+)?')
pat_end_block =   re.compile(r'^(#*\s*)?END (\w+)')
pat_comment = re.compile(r'^\s*#|^\s*$')
pat_property = re.compile(r'^(\w+?)\s*=\s*(.*)$')
pat_property_opt = re.compile('^ReplaceOption\[(\w+=\w+)\]=(.*)$')
pat_path_property  = re.compile(r'^(\S+?)::(\w+?)=(.*)$')

def isEndMarker(line, blockname):
    m = pat_end_block.match(line)
    return m and m.group(2) == blockname

def isComment(line):
    return pat_comment.match(line) is not None

def parse_old_plotfile(filename, hpath, section='PLOT'):
    """Parse a plot file of the old format and return a dict of the same format as the new plot files.

    Parameters
    ----------
    filename : str
        Filename of the .plot file. Should include the .plot extension.
    hpath : str
        The histogram path, usually with the format /AnalysisID/HistogramID.
    section : str
        Either PLOT or HISTOGRAM. Mostly just PLOT
    Returns
    -------
    plot_settings : dict
        Dict with of the same format and style as the new .plot files.

    Notes
    -----
    The SPECIAL and FUNCTION sections of a file cannot be parsed by the new matplotlib backend. 
    If these These will therefore be ignored.
    """
    if section not in ('PLOT', 'HISTOGRAM'):
        raise ValueError('Expected section to be PLOT or HISTOGRAM but got {}.'.format(section))

    # Decompose the histo path and remove the /REF prefix if necessary
    try:
        aop = AOPath(hpath)
    except:
        logging.debug("Found analysis object with non-standard path structure:", hpath, "... skipping")
        return {}

    ## Assemble the list of headers from any matching plotinfo paths and additional style files
    plotfile = aop.basepathparts()[0] + ".plot"

    pat_paths = {}
    ret = {}

    if not os.access(plotfile, os.R_OK):
        return
    startreading = False
    with open(plotfile) as f:
        msec = None
        for line in f:
            m = pat_begin_block.match(line)
            if m:
                tag, pathpat = m.group(2,3)
                # pathpat could be a regex
                if pathpat not in pat_paths:
                    try:
                        pat_paths[pathpat] = re.compile(pathpat)
                    except TypeError:
                        logging.debug("Error reading plot file for {}. Skipping.".format(plotfile))
                        return
                if tag == section:
                    m2 = pat_paths[pathpat].match(hpath)
                    if m2:
                        msec = m2
                        startreading = True
                        if section in ['SPECIAL']:
                            ret[section] = ''
                        continue
            if not startreading:
                continue
            if isEndMarker(line, section):
                startreading = False
                continue
            elif isComment(line):
                continue
            vm = pat_property.match(line)

            if vm:
                prop, value = vm.group(1,2)
                if msec:
                    oldval = value
                    try:
                        ## First escape backslashes *except* regex groups, then expand regex groups from path match
                        value = value.encode("string-escape")
                        value = re.sub("(\\\\)(\\d)", "\\2", value) #< r-strings actually made this harder, since the \) is still treated as an escape!
                        value = msec.expand(value)
                    except Exception as e: # TODO: bad exception handling
                        value = oldval #< roll back escapes if it goes wrong
                ret[section][prop] = texpand(value) #< expand TeX shorthands
            vm = pat_property_opt.match(line)
            if vm:
                prop, value = vm.group(1,2)
                ret[section]['ReplaceOption[' + prop + ']'] = texpand(value)

    return ret
