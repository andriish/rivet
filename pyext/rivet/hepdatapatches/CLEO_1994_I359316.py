import yoda
def patch(path, ao):
    # change bin center to CMS energy and add width
    if "CLEO_1994_I359316" in path and "d01" in path :
        for p in ao.points() :
            p.setXErrs(0.5)
            p.setX(10.58)
    return ao
