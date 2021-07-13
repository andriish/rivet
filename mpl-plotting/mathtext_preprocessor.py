# TODO: check that these macros have the correct replacement
macros = [['\GeV', r'\mathrm{Ge\!V}'], ['\TeV', r'\mathrm{Te\!V}'], 
          ['\pt', r'{\mathrm{p}_T}'], ['\pT', r'{\mathrm{p}_T}']]   # {} is added around to ensure that additional _ can be added in the label
def preprocess(s: str) -> str:
    """Convert convenient hepunits to mathtext.
    Preprocessor that converts certain commands in the input string. 
    Currently supported commands: `\GeV`, `\TeV`, `\pT`, `\pt`.
    This string will commonly be used for titles, x labels etc in figures in Rivet.
    Parameters
    ----------
    s : str
        Input string containing mathtext and one or multiple commands that will be converted.
    Returns
    -------
    s : str
        Output string that can be parsed and visualized by matplotlib. The string will be the same as s, except for the supported commands.
    """
    for macro, convert_to in macros:
        s = s.replace(macro, convert_to)
    
    return s