import yoda
def patch(path, ao):
    # make 3d scatters 2d
    if ("CMS_2017_I1485195" in path and type(ao) == yoda.core.Scatter3D) : 
        points = ao.points()
        newHist=yoda.Scatter2D()
        newHist.setPath(ao.path())
        for p in points:
            if "d06" in path :
                x     = p.y()
                xErrs = p.yErrs()
            else :
                x     = p.x()
                xErrs = p.xErrs()
            y     = p.z()
            newHist.addPoint(x,y)
            newHist.points()[-1].setXErrs(xErrs)
            for (key,value) in p.errMap().items() :
                newHist.points()[-1].setYErrs(value,key)
        ao=newHist
    return ao
