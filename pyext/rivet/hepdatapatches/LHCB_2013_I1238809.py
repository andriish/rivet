import yoda
def patch(path, ao):
    # bin widths
    if "LHCB_2013_I1238809" in path and "d01" in path:
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
