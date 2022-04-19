def patch(path, ao):
    # fix bin widths
    if "CRYSTAL_BALL_1990_I294492" in path :
        step=0.025
        if "d02" in path : step=0.05
        for p in ao.points() :
            p.setXErrs(step)
    return ao
