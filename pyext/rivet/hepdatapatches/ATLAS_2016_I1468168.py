
def patch(path, ao):
    if path == '/REF/ATLAS_2016_I1468168/d02-x01-y01':
      for p in ao.points:
          p.xErrs = (0.5, 0.5)
    return ao

