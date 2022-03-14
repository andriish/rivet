import yoda
def patch(path, ao):
    # missing bin widths (taken from fig 2 in paper)
    if "CDF_1997_I440101" in path :
        if "d02" in path or "d04" in path :
            bins = [5.,5.5,6.,6.5,7.,8.,9.,10.,12.,14.,17.,20.]
        elif "d03" in path or "d05" in path :
            bins = [5.,6.,7.,9.,12.,16.]
        else :
            return ao
        for (i,p) in enumerate(ao.points()) :
            p.setXErrs((p.x()-bins[i],bins[i+1]-p.x()))
    return ao
