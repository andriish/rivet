
def patch(path, ao):
    needs_patching = [ 
      '/REF/MARKII_1991_I295286/d01-x01-y01'
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setXErrs((0.1, 0.1))
    return ao

