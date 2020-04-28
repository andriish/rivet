
def patch(path, ao):
    needs_patching = [ 
      '/REF/UA5_1989_S1926373/d01-x01-y01',
      '/REF/UA5_1989_S1926373/d02-x01-y01',
      '/REF/UA5_1989_S1926373/d03-x01-y01',
      '/REF/UA5_1989_S1926373/d04-x01-y01',
      '/REF/UA5_1989_S1926373/d05-x01-y01',
      '/REF/UA5_1989_S1926373/d06-x01-y01',
      '/REF/UA5_1989_S1926373/d07-x01-y01',
      '/REF/UA5_1989_S1926373/d08-x01-y01',
      '/REF/UA5_1989_S1926373/d09-x01-y01',
      '/REF/UA5_1989_S1926373/d10-x01-y01',
      '/REF/UA5_1989_S1926373/d11-x01-y01',
      '/REF/UA5_1989_S1926373/d11-x01-y02',
      '/REF/UA5_1989_S1926373/d12-x01-y01',
      '/REF/UA5_1989_S1926373/d12-x01-y02',
    ]
    if path in needs_patching:
      bWidth = 0.5
      if 'd01' in path or 'd02' in path:
        bWidth = 0.1
      for p in ao.points():
          p.setErrs(1, (bWidth, bWidth))
    return ao

