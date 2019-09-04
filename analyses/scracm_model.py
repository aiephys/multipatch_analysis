"""
In sCRACM experiments, we photostimulate a large population of axon terminals while recording
postsynaptic responses. The intent is to be able to characterize the _strength_ of connectivity
between two cell populations, and in particular the relative strengths between different sets
of pre/postsynaptic populations. 

A few technical factors can get in the way of measuring these relative strengths accurately:
 - Differences in opsin expression between animals, between cells, and even between synapses
 - Differences in light exposure between synapses (because the photostimulation has nonuniform intensity)
 - Differences in sensitivity to photocurrent across synapses

The question we're exploring here is whether a stimuation protocol can be chosen that allows
us to work around the above limitations and yield more information about connection strength.

The model works as follows:
- A single postsynaptic cell is innervated by a large number of synapses
- Each synapse has a certain maximum synaptic strength
- Photostimulation causes activation of opsin channels based on an approximation of published
  intensity vs photocurrent relationships
- Each synapse receives a different fraction 0.5-1.0 of the photostimulation intensity to account
  for non-uniform illumination
- The level of opsin activation affects the probability of release via a nonlinear relationship.
  We don't really know what this relationship looks like (depends on calcium entry through opsin,
  binding to synaptic release mechanisms, etc.) so we model a few different possible relationships
  and explore the effect.
- Each synapse produces a postsynaptic response proportional to its release probability (rather
  than a stochastic/all-or-none release mechanism; assume this would all average out to the same
  result anyway)
- Synaptic responses are linearly summed to produce the final connection strength

"""

import numpy as np
from scipy.stats import norm, lognorm
import pyqtgraph as pg


# average strength of each type of connection to test
connection_types = [200e-6, 300e-6]

# what color to assign each connection type
type_colors = [(255, 0, 0, 100), (0, 255, 0, 100)]

# how many experiments to run, where each experiment yields one cell from 
# each connection type
n_experiments = 50

# how many presynaptic synapses connect to each cell
n_synapses = 100

# what stimulus intensities to test for each cell
power = 10**np.linspace(-2, 1, 100)

# distribution from which we select values that represent the
# level of opsin expression in each cell
#  - larger values of `scale` cause stronger opsin expression
#  - larger values of `s` cause more variance in opsin expression
per_cell_expression_distribution = lognorm(s=0.4, scale=4.0)


def sigmoid(x, xoffset, yoffset, slope, amp):
    return amp / (1.0 + np.exp(-slope * (x-xoffset))) + yoffset


class SCRACMModel:
    def __init__(self, syn_strength, n_synapses=100, chan_density=1.0):
        # average per-synapse strength
        self.syn_strength = syn_strength

        # number of synapses to model
        self.n_synapses = n_synapses

        # Light power is unevenly distributed across synapses
        self.stim_power_scale = np.random.uniform(0.5, 1.0, size=n_synapses)

        # fraction of open channels per synapse needed for 50% release probability
        self.release_threshold = lognorm.rvs(0.3, size=n_synapses) * 0.5 / chan_density

        # per-synapse strength scaling
        self.synapse_strength_scale = lognorm.rvs(0.5, size=n_synapses)

    def run(self, stim_power):
        """Return the sumulated evoked synaptic strength for a specific stimulus
        intensity.
        """
        # What fraction of channels are activated for each synapse
        channel_activation = self.channel_activation(stim_power * self.stim_power_scale)

        # Probability of vesicle release at each synapse
        release_prob = self.release_probability(channel_activation)

        # Total evoked synaptic strength
        total = (release_prob * self.synapse_strength_scale * self.syn_strength).sum()

        return total

    def channel_activation(self, stim_power):
        """Return the fraction of channels that open in response to stim_power (in mW/mm^2) at each synapse.
        """
        return 2**0.5 * sigmoid(stim_power, xoffset=0, yoffset=-0.5, amp=1, slope=1)**0.5

    def release_probability(self, chan_activation):
        """Return probability of releasing vesicle, given the number of open ion channels at each synapse.
        """
        return sigmoid(chan_activation / self.release_threshold, xoffset=0, yoffset=0, amp=1.0, slope=1)**16
        # return sigmoid(chan_activation / self.release_threshold, xoffset=0, yoffset=0, amp=1.0, slope=1)**16
        # return sigmoid(chan_activation, xoffset=self.release_threshold, yoffset=0, amp=1.0, slope=30)



plt = pg.plot(
    labels={'left': ('evoked synaptic strength', 'V'), 'bottom': 'stimulus intensity (mW/mm^2)'}, 
    title="Dose-response curves for %d connection types, each with %d cells" % (len(connection_types), n_experiments),
)
plt.setLogMode(x=True)

# Run all experiments, taking note of final connection strength and slope for each
max_strength = [[], []]
end_slope = [[], []]
for i in range(n_experiments):
    for j, (strength, color) in enumerate(zip(connection_types, type_colors)):
        # for each cell to test, we pick a random value to represent the variance in
        # opsin expression between cells and between animals.
        chan_expression = per_cell_expression_distribution.rvs()

        # initialize model
        model = SCRACMModel(syn_strength=strength, n_synapses=n_synapses, chan_density=chan_expression)

        # run model once for each stim power
        tot = [model.run(stim_power=p) for p in power]

        # record some results
        max_strength[j].append(tot[-1])
        end_slope[j].append(tot[-1]-tot[-2])

        # plot the dose-response curve for this cell
        plt.plot(power, tot, pen=color)


# beeswarm plot of maximum connection strength for each cell
#   (seems to be the best way to determine relative connection strength)
plt = pg.plot(
    title="Response strengths at maximum stimulus intensity",
    labels={'left': ('evoked synaptic strength', 'V'), 'bottom': 'connection type'}
)
plt.getAxis('bottom').setTicks([[(i, 'type %d'%i) for i in range(len(connection_types))]])
for i,color in enumerate(type_colors):
    y = np.array(max_strength[i])
    x = pg.pseudoScatter(y, bidir=True)
    x = 0.2 * x / x.max()
    plt.plot(x+i, y, pen=None, symbol='o', symbolBrush=color, symbolPen=None)
x = np.array([0, 1])
y = np.array(max_strength).mean(axis=1)
h = np.array(max_strength).std(axis=1)
err = pg.ErrorBarItem(x=x, y=y, height=h, width=0.1)
plt.addItem(err)


# Plot the slope at the end of each dose-response curve versus max connection strength
# 
plt = pg.plot(labels={'left': 'final normalized dose-response slope', 'bottom': 'evoked synaptic strength'})
for i,color in enumerate(type_colors):
    plt.plot(max_strength[i], np.array(end_slope[i]) / np.array(max_strength[i]), pen=None, symbol='o', symbolBrush=color, symbolPen=None)
plt.setLogMode(y=True)


# Plot the relationship between stimulus intensity and opsin channel activation
open_ch = model.channel_activation(power)
plt = pg.plot(power, open_ch, labels={'bottom': 'stim power (mW/mm^2)', 'left': 'opsin channel activation'})
plt.setLogMode(x=True)


# For the last model that was generated, plot the relationships between opsin channel activation 
# and release probability for all synapses in the model
chan_act = np.linspace(0, 1, 100)
plt = pg.plot(labels={'left': 'release probability', 'bottom': ' opsin channel activation'})
rprob = model.release_probability(chan_act[:, None])
for i in range(rprob.shape[1]):
    plt.plot(chan_act, rprob[:,i])

