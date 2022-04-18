def patch(path, ao):
    # fix missing bin widths
    if "TASSO_1983_I191417" in path and ("d01" in path or "d02" in path):
        for p in ao.points() :
            p.setXErrs(0.05)
    return ao
