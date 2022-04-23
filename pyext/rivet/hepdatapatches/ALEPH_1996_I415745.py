import yoda
def patch(path, ao):
    # remove average etc bins and make into separate dists
    if "ALEPH_1996_I415745" in path and ("d01" in path or "d02" in path):
        if "d01" in path :
            bins=[.1,0.15,.2,.3,.4,1.]
        else :
            bins=[.3,.6,.9,1.2,1.5]
        output=[]
        for i, p in enumerate(ao.points()):
            if i < len(bins)-1:
                x  = 0.5*(bins[i+1]+bins[i])
                dx = 0.5*(bins[i+1]-bins[i])
                p.setX(x)
                p.setXErrs(dx)
            else :
                newHist=yoda.Scatter2D()
                newHist.setPath(ao.path().replace("x01","x0%s"%(i-len(bins)+3)))
                if "d01" in path :
                    p.setX(0.65)
                    p.setXErrs(0.35)
                elif i==4 :
                    p.setX(5.75)
                    p.setXErrs(4.25)
                elif i==5 :
                    p.setX(5.15)
                    p.setXErrs(4.85)
                elif i==6 :
                    p.setX(5.3)
                    p.setXErrs(4.7)
                newHist.addPoint(p)
                output.append(newHist)
        imax=len(ao.points())
        for i in range(len(bins)-1,imax) :
            ao.rmPoint(i)
        output.append(ao)
        return output
    return ao
