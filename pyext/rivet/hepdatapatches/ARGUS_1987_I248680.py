def patch(path, ao):
    # fix missing bin widths
    if "ARGUS_1987_I248680" in path:
        if "d01" in path :
            step=0.05
        elif "d02" in path :
            step=0.15
        elif "d03" in path :
            step=0.2
        elif "d04" in path :
            step=0.4
        else :
            return ao
        for p in ao.points() :
            p.setXErrs(step)
    return ao
