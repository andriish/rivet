def patch(path, ao):
    # fix dodgy bin widths
    if "BELLE_2009_I815978" in path and "d31" in path:
        for p in ao.points() :
            if p.x()<1.8 :
                p.setXErrs(0.01)
    return ao
