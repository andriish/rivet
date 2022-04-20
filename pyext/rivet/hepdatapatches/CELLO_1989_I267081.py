def patch(path, ao):
    # fix missing bin widths
    if "CELLO_1989_I267081" in path:
        for p in ao.points() :
            p.setXErrs(0.15)
    return ao
