def patch(path, ao):
    # fix missing bin widths
    if "TPC_1987_I246557" in path and "d01" in path:
        step=0.05
        for p in ao.points() :
            if p.x()==2.5 : step=0.1
            p.setXErrs(step)
            if p.x()==2.2 : p.setX(2.25)
    return ao
