import yoda
def patch(path, ao):
    # add bin widths
    if ( "ALICE_2021_I1898832" in path and
         ("d01" in path or "d02" in path or "d03" in path) ) : 
        points = ao.points()
        if "d01" in path :
            edges=[0.,0.5,1.,2.,3.,4.,5.,7.,10.,15.,20.,25.,30.,40.]
        elif "d02" in path:
            edges=[-0.9,0.9]
        elif "d03" in path:
            edges=[0.,0.2,0.5,0.9]
        for i, p in enumerate(ao.points()):
            p.setXErrs((p.x()-edges[i],edges[i+1]-p.x()))
    return ao
