import yoda
def patch(path, ao):
    # dists integrated over bin, undo this
    if "CMS_2013_I1225274" in path :
        ihist = int(path.split("/")[-1].split("-")[0][1:])
        # hists we don't need to touch
        if ihist == 1 :
            return ao
        factor=1.
        # divide by rapidity bining for double differential and undo folding +/- rapidity
        if(ihist>=3 and ihist<=22) :
            factor *= 0.5/0.4
        # undo folding +/- rapidity
        elif(ihist>=23 and ihist<=25) :
            factor *= 0.5
        for p in ao.points():
            xerrs = p.xErrs()
            width= xerrs[0]+xerrs[1]
            p.setY(p.y()*factor/width)
            for (key,value) in p.errMap().items() :
                temp=(value[0]*factor/width,value[1]*factor/width)
                p.setYErrs(temp,key)
    return ao
