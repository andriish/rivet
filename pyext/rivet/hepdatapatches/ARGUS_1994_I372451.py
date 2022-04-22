def patch(path, ao):
    # fix bin value and width due limit bin in hepdata
    if "ARGUS_1994_I372451" in path and "d02" in path:
        for p in ao.points() :
            p.setX(2.1)
            p.setXErrs(0.2)
    return ao
