def patch(path, ao):
    # fix missing bin widths
    if "ARGUS_1989_I267759" in path and "d02" not in path and "d07" not in path:
        step=0.1
        if "d01" in path : step=0.05
        for p in ao.points() :
            p.setXErrs(step)
    return ao
