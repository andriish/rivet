"""This module creates a rivet-style plot of d10-x01-y01 using mpl."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import data  # Import xlow, xhigh, val, errminus, errplus of data
import mc1  # Import xlow, xhigh, val, errminus, errplus of mc1
from style import style  # Import the style dictionary used for rcparams

plt.style.use(style)
fig, axes = plt.subplots(2, 1, sharex=True,
                         gridspec_kw={'height_ratios': [2, 1]})

# Set axis limits
xmin = min(data.xlow)
xmax = 50

y1min = 6e-7
y1max = 3
y2min = 0.5
y2max = 1.4999999  # Avoids ytick label at 1.5

axes[0].set_xlim(xmin, xmax)
axes[0].set_ylim(y1min, y1max)
axes[1].set_ylim(y2min, y2max)
axes[0].set_xscale('log')
axes[0].set_yscale('log')

# Set axis ticks
x_major = mpl.ticker.LogLocator(base=10.0)
axes[0].xaxis.set_major_locator(x_major)
y_major = mpl.ticker.LogLocator(base=10.0)
axes[0].yaxis.set_major_locator(y_major)

# Top plot
# Data line
axes[0].hlines(data.val, data.xlow, data.xhigh, 'k')
x_data = (data.xlow + data.xhigh)/2
y_data = data.val
axes[0].plot(x_data, y_data, 'ko')
axes[0].vlines(x_data, (data.val - data.errminus),
               (data.val+data.errplus), 'k', zorder=3)
# mc1 line
x_mc1 = np.append(mc1.xlow, mc1.xhigh[-1])
y_mc1 = np.insert(mc1.val, 0, mc1.val[0])
axes[0].plot(x_mc1, y_mc1, 'r', drawstyle='steps-pre', solid_joinstyle='miter')

# Bottom plot
x_ratio = x_mc1
y_ratio = y_mc1/np.insert(y_data, 0, y_data[0])
axes[1].plot(x_ratio, y_ratio, 'r', drawstyle='steps-pre', zorder=1,
             solid_joinstyle='miter')
axes[1].hlines(1, xmin, xmax, 'k', zorder=2)
axes[1].plot(x_data, np.ones(len(x_data)), 'ko', zorder=3)
axes[1].vlines(x_data, (data.val - data.errminus)/data.val,
               (data.val+data.errplus)/data.val, 'k', zorder=3)
axes[1].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))

# Axes labels
axes[1].set_xlabel('$p_\perp$ [GeV]')
axes[0].set_ylabel('$1/N_\mathrm{ev} \, $1$/2\pi{}p_\perp \, \mathrm{d}\sigma/\mathrm{d}\eta\mathrm{d}p_\perp$',
                   loc='top')
axes[1].set_ylabel('MC/Data')
axes[0].set_title('Charged particle $p_\perp$ at 7 TeV, track $p_\perp > 500$ MeV, for $N_\mathrm{ch} \geq 1$',
                  loc='left')

# Legend
# TODO: Find a better way to implement the custom Data legend graphic


class AnyObject(object):
    """Necessary for custom legend handler."""
    pass


class AnyObjectHandler(object):
    """Creates custom legend handler for Data."""

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        patch = mpl.patches.Circle(
            (x0+width/2, y0+height/2), (2.5)**0.5, facecolor='black')
        handlebox.add_artist(patch)
        patch = mpl.patches.Rectangle(
            (-0.4, 3), width+0.8, 0.8, facecolor='black')
        handlebox.add_artist(patch)
        patch = mpl.patches.Rectangle(
            (width/2-0.4, 0), 0.8, height, facecolor='black')
        handlebox.add_artist(patch)
        return patch


handler_mc1 = mpl.lines.Line2D([], [], color='red')
axes[0].legend([AnyObject(), handler_mc1], ['Data', 'mc1'],
               handler_map={AnyObject: AnyObjectHandler()}, loc=[0.52, 0.77])

plt.savefig('mpl_d10-x01-y01.pdf', dpi=500)
