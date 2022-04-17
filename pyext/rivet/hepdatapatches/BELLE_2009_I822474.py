def patch(path, ao):
    # fix dodgy bin
    if "BELLE_2009_I822474" in path and "d01" in path:
        for p in ao.points() :
            if p.x()==0.99 :
                p.setXErrs(0.01)
    return ao
