import yoda
def patch(path, ao):
    # divide by bin width to make a differential distribution
    if "L3_1998_I447945" in path :
        for p in ao.points() :
            width=p.xErrs().minus+p.xErrs().plus
            p.setY(p.y()/width)
            for (key,value) in p.errMap().items() :
                p.setYErrs((value[0]/width,value[1]/width),key)
    return ao
