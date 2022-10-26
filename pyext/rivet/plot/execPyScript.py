import os

def execute(pyScript):
  if not os.path.isfile(pyScript):
    raise IOError("Python script not found!")

  # supply globals() argument to add imports to global scope
  return exec(open(pyScript).read(), globals())
  
