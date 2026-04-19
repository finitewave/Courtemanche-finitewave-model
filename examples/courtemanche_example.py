"""
Example script to run a 0D model simulation and plot the results.

This script sets up a simple stimulation protocol, runs the simulation,
and plots the membrane potential over time.
"""

import numpy as np
import matplotlib.pyplot as plt

from implementation.courtemanche_0d import Courtemanche0D, Stimulation


stimulations = []
stimulations.append(Stimulation(t_start=100, duration=2, amplitude=20.0))

t_max = 600  # Total simulation time in ms

model = Courtemanche0D(dt=0.01, stimulations=stimulations)
model.run(t_max=t_max)

time = np.arange(0, t_max, model.dt)

fig = plt.figure()
plt.plot(time, model.history['u'], lw=2)
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (u)')
plt.title('0D Courtemanche Simulation')
plt.grid(which='major')
plt.show()

# fig.savefig("courtemanche_ap.png", dpi=300)
