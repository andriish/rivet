import yoda
def patch(path, ao):
    # fix hist, really a 2D histo with average of x values added
    if ("CMS_2015_I1342266" in path and
        ("d01" in path or "d02" in path or "d03" in path)) : 
        points = ao.points()
        newHist=yoda.Scatter2D()
        newHist.setPath(ao.path())
        for p in points:
            x     = p.x()
            xErrs = p.xErrs()
            y     = p.z()
            yErrs = p.zErrs()
            newHist.addPoint(x,y,xErrs,yErrs)    
        ao=newHist
    return ao
