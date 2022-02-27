import yoda
def patch(path, ao):
    # add dummy bin width for plotting
    if ("CMS_2012_I1113442" in path and "d01" in path) :
        for p in ao.points():
            p.setXErrs((0.5,0.5))
    return ao
