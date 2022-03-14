import yoda
def patch(path, ao):
    # missing bin widths
    if "CDF_1997_I440446" in path :
        if "d01" in path:
            bins = [1799.5,1800.5]
        else :
            bins = [5.,5.5,6.,6.5,7.,8.,9.,10.,12.,14.,17.,20.]
        for (i,p) in enumerate(ao.points()) :
            p.setXErrs((p.x()-bins[i],bins[i+1]-p.x()))
    return ao
