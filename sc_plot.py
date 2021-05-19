import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import os

tf = 1e6
timeline = np.linspace(0, tf, 10000)
G = 6.674e-11
txt_nmass = np.genfromtxt('Nmass.txt', delimiter=',')
mstar = txt_nmass[0][8]
mu_paramater = G*mstar

orbital_elements = [orb_ele for orb_ele in os.listdir('.') if orb_ele.startswith('OrbitalElements')]

temp_energy = np.array([])

for energy_index in orbital_elements:
    txt_orb_elements = np.genfromtxt(energy_index, delimiter=',')
    semimajor = txt_orb_elements[:,2]

#    print(mu_paramater.shape,semimajor.shape)
    specific_orbital_energy = -mu_paramater/(2*semimajor)
    temp_energy = np.append(temp_energy, np.sum(specific_orbital_energy))

energy_ini = temp_energy[0]
plt.plot(timeline, (temp_energy - energy_ini)/energy_ini)
plt.savefig('../plot.png')
