# Required packages
### rivet_makeyaml
- ruamel.yaml (might be replaced by pyyaml in the future)
### testing
- pytest
- pytest-mpl (currently not used but mihgt be added)
### 2D histogram
- matplotlib
- numpy (a dependency of matplotlib)

# Mathtext-compatible .plot files

## Fixed issues
- Changed `\text` to `\mathrm`.
- Added curly braces to spots where it is missing from `\frac`, `\sqrt`, `\text`.
- Changed `\ge` to `\geq`, `\le` to `\leq`.
- Fixed edge cases where converting `\text` to `\mathrm` does not work.
- Remove `^` from `\perp^` in `analyses/pluginCMS/CMS_2014_I1305624.plot`.
- Replace `\unit` and `\si` with alternatives, e.g., `\unit{10}{\GeV}` is converted to `10~\GeV`.
- Spaces inside former `\text` commands have been replaced by a `~`, since `\mathrm` ignores spaces while `\text` does not.
- Enclose LaTeX commands such as `\%` with $ (i.e. math mode). Otherwise, mathtext will interpret `\` as plain text. 
- Replace deprecated `\rm` command with `\mathrm`
- Remove superfluous `}` in some labels that caused LaTeX syntax errors
- Replace `\PQt` with `\mathrm{t}`, since `\PQt` is a non-standard command in LaTeX and is not recognized by matplotlib's mathtext.

## Remaining issues
- `\micro` in LaTeX will give the same output as `$\mathrm{\mu}$` in mathtext. However, `$\mathrm{\mu}$` in LaTeX will give a [non-sensical output](https://tex.stackexchange.com/questions/569676/applying-mathrm-on-mu-leads-to-strange-symbol). This can therefore not be changed until the backend has been switched completely to mathtext.
- `\large` has no equivalent command in mathtext. Instead, one must pass a `fontsize` parameter to matplotlib. The details for how such a parameter should be passed to matplotlib has not been decided yet. As such, this cannot be changed for now.
- `madgraph_a` cannot be parsed by LaTeX but can be parsed by mathtext and will be written as plaintext. Is this intended?

## Commands that will not be converted
`\GeV`, `\TeV`, `\pt` and `\pT`, and potentially some other commonly used simple commands will not be converted to something else in this merge request. The idea is instead to preprocess each label before it is passed to matplotlib and replace these commands with something that can be parsed by mathtext.

# Installing fonts
## How to make matplotlib find the relevant fonts 
For Rivet, `URW Palladio L` and `PazoMath` are the most relevant ones but this method works for any installed font.
1. Install texlive-fonts-recommended to get many fonts (not required)
2. Font files will be in /usr/share/. The font path can probably be found by executing `find | grep fontname` inside /usr/share/. For pazomath: `find | grep mathpazo`. 
    - pazomath should be located in `/usr/share/texlive/texmf-dist/fonts/type1/public/mathpazo`.
3. Convert the pfb files into ttf files 
    - Might not be necessary if matplotlib can read pfb files. However, the method for this is unknown. Use e.g. https://convertio.co/ (not sure if website is safe).
4. Make matplotlib find the ttf files. This can be done in many ways:
    - Copy the ttf files into a place where matplotlib looks for fonts, e.g. by creating the directory /usr/share/fonts/truetype/fontname (if they are ttf files. If matplotlib can read pfb files, move them to /usr/share/fonts/type1, including afm files). 
    - The method below is not tested but would be a better option, since the user will not need to move font files etc.
```python
import matplotlib as mpl
import matplotlib.font_manager as fm

fe = fm.FontEntry(fname='your custom ttf file path', name='your custom ttf font name')
fm.fontManager.ttflist.insert(0, fe) # or append is fine
mpl.rcParams['font.family'] = 'your custom ttf font name'
```
5. Clear the font cache by removing `~/.cache/matplotlib/` and running `fc-cache -fv` (at least on Ubuntu). 
6. From matplotlib, change the font of mathtext and/or the regular text font. As an example:
    - Change regular font to `URW Palladio L`: TODO: double check this command.

```python
mpl.rcParams.update({'font.family': 'sans-serif', 'font.serif': 'URW Palladio L'})
```

- Change mathtext font to `PazoMath`:

```python
mpl.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': 'PazoMath', 'mathtext.bf': 'PazoMath:bold', 'mathtext.it': 'PazoMath:italic'})
```
 - I still do not know how to get calligraphy letters with PazoMath and how to get old style fonts. 

# New format for .plot and .dat files
The current format of .plot and .dat files has been replaced by a YAML-based format which makes it easier to read and parse.

The content will be replaced to better fit the new matplotlib backend. It currently contains these sections (names not final):
- **style**: A predefined style, such as `default`, which will display the plot in a certain way by using a specific font etc. Mainly (or only) controlled by rcparams. 
- **rcparams**: matplotlib rcparams that will change the look of the figure. These will be applied after the style rcparams are applied.
- **plot features**: contains all settings that were previously allowed, such as Title, YLabel etc.
- **histograms (only in .dat files)**: a dictionary of histograms with their names, values and other properties, including a large string corresponding to the histogram written in the YODA format.

- New arguments and removed arguments in files, such as 
    - `SPECIAL` section (removed)
    - `FUNCTION` section (removed)
    - `2DType` in plot features (name not final, added)
    - `2DIndividual` in plot features (name not final, added)
    - `2DRatioColormap` in plot features (name not final, added)
    - `2DColormap` in plot features (name not final, added)
- The format of .dat and .plot files (including examples)
- Difference between .plot and .dat files (histograms do not exist in .plot files and .dat files they have different use cases)

# Updated rivet-cmphistos
`rivet-cmphistos` has now been replaced by `rivet_makeyaml` (name not final, might be renamed to `rivet_mkdat`).

## New behavior for 2D histograms
Due to limitations to the old plotting backend (i.e. `make-plots`), multiple 2D histograms could not be plotted in the same figure. To circumvent this, `rivet-cmphistos` would output 2D histograms with the same path (i.e. `analysis ID/histogram ID`) into separate .dat files. 
The new plotting backend can plot multiple 2D histograms in the same figure (see the [2D histograms section](#2D-histograms) below). No hacky behavior is therefore needed inside `rivet_makeyaml` and all 2D histograms will therefore be included in the same output yaml file, just like how 1D histograms behave.

## Updated flags for rivet-cmphistos
- `rivet-cmphistos` has been renamed `rivet_makeyaml` (name not final) and will mainly be called as a python function by `rivet-mkhtml`. In other words, it might not remain as a stand-alone bash command that can be called. However, if this is desired, it is easy to add a thin wrapper around rivet_makeyaml to make a CLI.
- Since everything will be controlled by `rivet-mkhtml`, all command line arguments will be passed to `rivet-mkhtml`. A subset of these will then be passed to the `make_yamlfiles` function.
- See [this table](https://docs.google.com/spreadsheets/d/1GUpjXIZToN0vr4dkPfJBvTReHq9j3AqRx5v5Yeel95U/edit?usp=sharing) of all new and old arguments. TODO move this table to this document as a md table. 

## YAML file parsing (only important for developers)
All IO of yaml files is in yamlio.py. This way, if we need to change backend, only that file and the import needs to be changed there.

# New plot 
# 2D histograms
There are multiple ways for how 2D histograms can be plotted:

## Output format
In this section, "output file" means an image, i.e., a png, pdf, svg etc.

### Output files
- Save each 2D histogram and ratio plot to a separate output file (default)
- Save all 2D histograms and ratio plots in the same output file as a grid of subplots

### Visualization mode
- Plot 2D histograms and ratio plots as heatmaps on the x-y plane, with the color of each bin indicating its "height" (default)
- Plot 2D histograms and ratio plots as surface plots.

The 2 output file modes can be combined with the 2 visualization modes. 

## Code location
The code for plotting 2D histograms is in `rivet_plot2d.py` (name not final). The functions inside this file will be called by `rivet_plot.py`.

Parts of the code will be moved to yoda and become a part of its plotting API. these will then become functions that will likely take yoda histograms and plot unformatted plots. Rivet then builds upon this by using the yoda plotting API and adding formatting options using both the plotting API and some raw matplotlib. 

## Examples
```bash
rivet-mkhtml mc1.yoda mc2.yoda
```
