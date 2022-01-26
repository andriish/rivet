import yoda
def patch(path, ao):
    # fix hist, this one is a mess really 3 different 2d histos
    if ("LHCB_2013_I1244315" in path and ("d01" in path or "d02" in path)) :
        points = ao.points()
        newHists=[]
        for i in range(1,4) :
            newHists.append(yoda.Scatter2D())
            newHists[-1].setPath(ao.path().replace("x01","x0%s"%i))
        for p in points:
            x     = p.x()
            xErrs = p.xErrs()
            y     = p.z()
            yErrs = p.zErrs()
            iloc = int(p.y()-1)%3
            newHists[iloc].addPoint(x,y,xErrs,yErrs)    
        return newHists
    else :
        return ao
