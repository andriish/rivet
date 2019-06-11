
def patch(path, yodaobject):
    if path == '/REF/ATLAS_2016_I1468168/d02-x01-y01':
      for p in yodaobject.points:
          p.xErrs = (0.5, 0.5)
    return yodaobject

