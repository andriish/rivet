import yoda
def patch(path, ao):
    # dists integrated over bin, undo this
    if "LHCB_2014_I1280929" in path :
        ihist = int(path.split("/")[-1].split("-")[0][1:])
        if ihist >= 3 and ihist <=8 :
            for p in ao.points():
                xerrs = p.xErrs()
                width= xerrs[0]+xerrs[1]
                p.setY(p.y()/width)
                for (key,value) in p.errMap().items() :
                    temp=(value[0]/width,value[1]/width)
                    p.setYErrs(temp,key)
    return ao
