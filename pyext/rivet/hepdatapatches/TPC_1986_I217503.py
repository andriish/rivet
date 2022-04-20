def patch(path, ao):
    # fix missing bin widths
    if "TPC_1986_I217503" in path:
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
