def patch(path, ao):
    # fix bin widths
    if ("ATLAS_2018_I1622737" in path and
        ("d01" in path or "d02" in path or "d03" in path or "d04" in path or
         "d05" in path or "d06" in path or "d07" in path)) :
        if "d01" in path or "d03" in path :
            bins=[8.,9.,10.,11.,12.,13.,14.,16.,18.,22.,30.,40.]
        elif "d02" in path or "d04" in path :
            bins=[8.,10.,12.,16.,22.,40.]
        elif "d05" in path :
            bins=[0.,2.,4.,6.,8.,10.,14.,20.,40.]
        else :
            bins=[0.,6.,10.,15.,40.]
        for i in range(0,len(ao.points())) :
            x = ao.points()[i].x()
            ao.points()[i].setXErrs((x-bins[i],bins[i+1]-x))
    return ao
