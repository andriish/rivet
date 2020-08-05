# bin widths
def patch(path, ao):
    if "L3_2004_I652683" in path and "d59" in path:
        for p in ao.points():
            p.setXErrs(1)
    return ao
