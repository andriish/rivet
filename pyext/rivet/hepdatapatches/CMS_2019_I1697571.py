def patch(path, ao):
    # fix bin widths
    if "CMS_2019_I1697571" in path:
        if "d01" in path :
            bins =[7,15,20,50]
        else :
            bins = [7,15,50]
        for i,p in enumerate(ao.points()):
            p.setXErrs((p.x()-bins[i],bins[i+1]-p.x()))
    return ao
