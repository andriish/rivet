#!/usr/bin/env python
"""%prog [OPTIONS] <AIDAFILE> [...]

Strip specified bins from data sets. Histgrams not specified will be passed
through without any chopping. Bins to be kept can be specified on command
line via `-b' options. The format is
    -b AIDAPATH:start:stop
where start and stop are x values contained in the first and last bins,
respectively, that should be kept. They need not to be the center but must
only lie somewhere in the bin's x-range.

Example:
    %prog -b /ALEPH_1996_S3486095/d03-x01-y01:0.095:0.27 out.aida
This will give you the all bins of the ALEPH 1-T distribution that are
between the bins that contain the x-values 0.095 and 0.27 .
"""

import os
import sys
import logging

import lighthisto
## Make "sorted" a builtin function on Python < 2.4
if not 'sorted' in dir(__builtins__):
    def sorted(iterable, cmp=None, key=None, reverse=None):
        rtn = iterable
        rtn.sort(cmp)
        return rtn

## Add logging.log if needed
if not 'log' in dir(logging):
    def _logit(level, msg):
        l = logging.getLogger()
        l.log(level, msg)
    logging.log = _logit


## Try to load faster but non-standard cElementTree module
try:
    import xml.etree.cElementTree as ET
except ImportError:
    try:
        import cElementTree as ET
    except ImportError:
        try:
            import xml.etree.ElementTree as ET
        except:
            sys.stderr.write("Can't load the ElementTree XML parser: please install it!\n")
            sys.exit(1)


if __name__ == "__main__":
    from optparse import OptionParser, OptionGroup
    parser = OptionParser(usage=__doc__)

    parser.add_option("-b", "--bins",
                      action="append",
                      help="Specify a histogram and bin range that is to be"
                           " kept. The format is `AIDAPATH:start:stop'.")
    parser.add_option("-o", "--out",
                      dest="outdir",
                      help="output directory (default: %default)")

    verbgroup = OptionGroup(parser, "Verbosity control")
    verbgroup.add_option("-V", "--verbose", action="store_const",
                         const=logging.DEBUG, dest="LOGLEVEL",
                         help="print debug (very verbose) messages")
    verbgroup.add_option("-Q", "--quiet", action="store_const",
                         const=logging.WARNING, dest="LOGLEVEL",
                         help="be very quiet")
    parser.set_defaults(bins=[],
            outdir=".",
            LOGLEVEL=logging.INFO)
    opts, args = parser.parse_args()

    ## Configure logging
    try:
        logging.basicConfig(level=opts.LOGLEVEL,
                format="%(levelname)s: %(message)s")
    except:
        logging.getLogger().setLevel(opts.LOGLEVEL)
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        logging.getLogger().addHandler(h)

    if len(args) == 0:
        sys.stderr.write("Must specify at least one AIDA histogram file\n")
        sys.exit(1)

    bindefs = {}
    for bd in opts.bins:
        try:
            path, low, high = bd.split(":")
        except:
            sys.stderr.write("Problem parsing bin definition `%s'" % (bd))
            sys.exit(1)
        if low == "":
            low = None
        else:
            low = float(low)
        if high == "":
            high = None
        else:
            high = float(high)
        bindefs[path] = (low, high)

    for aidafile in args:
        if not os.access(aidafile, os.R_OK):
            logging.error("%s can not be read" % aidafile)
            break

        base, ext = os.path.splitext(os.path.basename(aidafile))
        chopfile = os.path.join(opts.outdir, base + "-chop" + ext)
        outhistos = []

        tree = ET.parse(aidafile)
        for dps in tree.findall("dataPointSet"):
            thishist = lighthisto.Histo.fromDPS(dps)
            if thishist.fullPath() in bindefs.keys():
                outhistos.append(thishist.chop(bindefs[thishist.fullPath()]))
            else:
                outhistos.append(thishist)
        out = open(chopfile, "w")
        out.write('<?xml version="1.0" encoding="ISO-8859-1" ?>\n')
        out.write('<!DOCTYPE aida SYSTEM "http://aida.freehep.org/schemas/3.3/aida.dtd">\n')
        out.write('<aida version="3.3">\n')
        out.write('  <implementation version="1.1" package="FreeHEP"/>\n')
        out.write("\n\n".join([h.asAIDA() for h in sorted(outhistos)]) + "\n")
        out.write("</aida>\n")
        out.close()
