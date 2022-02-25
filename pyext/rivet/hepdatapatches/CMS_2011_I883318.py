import yoda
def patch(path, ao):
    # fix hist, bin index rather than limits used
    if "CMS_2011_I883318" in path and "d02-x01-y01" in path:
        bins=[5.,10.,13.,17.,24.,30.]
        for (i,p) in enumerate(ao.points()) :
            x  = 0.5*(bins[i+1]+bins[i])
            dx = 0.5*(bins[i+1]-bins[i])
            p.setX(x)
            p.setXErrs((dx,dx))
    return ao
