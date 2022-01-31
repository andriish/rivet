import yoda
def patch(path, ao):
    # fix hist, really a 2D histo with average of x values added
    if ("CMS_2011_I878118" in path and
        ("d05" in path or "d06" in path or "d07" in path)) : 
        points = ao.points()
        newHist=yoda.Scatter2D()
        newHist.setPath(ao.path())
        for p in points:
            x     = p.x()
            xErrs = p.xErrs()
            y     = p.z()
            newHist.addPoint(x,y)
            newHist.points()[-1].setXErrs(xErrs)
            for (key,value) in p.errMap().items() :
                newHist.points()[-1].setYErrs(value,key)
        ao=newHist
    return ao
