def patch(path, ao):
    # fix missing bin widths
    if "ARGUS_1988_I262713" in path:
        if "d04" in path :
            step=0.2
        else :
            step=0.1
        for p in ao.points() :
            if p.x()>2.75 :
                p.setXErrs(0.2)
            elif "d04" in path and p.x()>2.65 :
                p.setXErrs(0.3)
            else :
                p.setXErrs(step)
    return ao
