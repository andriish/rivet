import os
from multiprocessing import Pool

def execute(pyScript):
  if not os.path.isfile(pyScript):
    raise FileNotFoundError("Python script not found!")

  # supply globals() argument to add imports to global scope
  exec(open(pyScript).read(), globals())
  

def executePyScripts(pyScripts): 
  from multiprocessing import Pool
  p = Pool(processes=4)
  p.map(execute, pyScripts)
  
