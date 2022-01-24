def patch(path, ao):
    # fix bin widths
    if "ALICE_2019_I1735351" in path and "d02" in path :
        bins = [0.,1.,2.,3.,4.,5.,7.,10.]
        for i in range(0,len(ao.points())) :
            x = ao.points()[i].x()
            ao.points()[i].setXErrs((x-bins[i],bins[i+1]-x))
    return ao
