import yoda_plot1d as yp
import yoda
import matplotlib.pyplot as plt
import numpy as np

# Create three Histo1D objects
rng = np.random.default_rng(seed=1)
nevents = 10_000
h1, h2, h3 = yoda.Histo1D(10, -1, 5), yoda.Histo1D(10, -1, 5), yoda.Histo1D(10, -1, 5)
for h in (h1, h2, h3):
    for event in rng.normal(1.5, 1, size=(nevents, 2)):
        h.fill(event[0], event[1])
h1, h2, h3 = h1.mkScatter(), h2.mkScatter(), h3.mkScatter()  # TODO: Move into yoda 

# Plot objects using "yoda"
ax = yp.plot_hist((h1, h2, h3), ErrorBars=True, use_reference=True)  # TODO: Change to this
plt.savefig('yoda1D_hist.jpg')

ax = yp.plot_ratio((h1, h2, h3), ErrorBars=True)
plt.savefig('yoda1D_ratio.jpg')
