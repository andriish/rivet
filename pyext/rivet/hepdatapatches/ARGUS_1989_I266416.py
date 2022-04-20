def patch(path, ao):
    # fix missing bin widths
    if "ARGUS_1989_I266416" in path:
        for p in ao.points() :
            if "d01" in path :
                p.setXErrs(0.05)
            else :
                if p.x()>1.2 and p.x()<2.4 :
                    p.setXErrs(.1)
                elif p.x()<1. :
                    p.setXErrs((0.15,0.25))
                elif p.x()<3. :
                    p.setXErrs(0.2)
                else :
                    p.setXErrs(0.35)
    return ao
