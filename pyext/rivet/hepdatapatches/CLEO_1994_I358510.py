def patch(path, ao):
    # fix missing bin widths
    if "CLEO_1994_I358510" in path and "d01" not in path:
        for p in ao.points() :
            p.setXErrs(0.05)
    return ao
