def patch(path, ao):
    if "MARKIII_1987_I244856" in path:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
