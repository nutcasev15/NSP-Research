# OpenMC 0.13.3 Python Post Processing Script for Reflector BOL Keff Study


############### Library Imports
import pickle
import itertools
import matplotlib.markers as mrk
import matplotlib.pyplot as plt
import matplotlib.ticker as tck


############### Load Results from Pickle File
Pikname = 'BOL-Ref-0'
with open(Pikname + '.pckl', 'rb') as f:
    Rings = pickle.load(f)
    pow = pickle.load(f)
    Ref_thic = pickle.load(f)
    MF_Mod = pickle.load(f)
    Mname = pickle.load(f)
    Rname = pickle.load(f)
    core_rad = pickle.load(f)
    core_hei = pickle.load(f)
    ref_mass = pickle.load(f)
    core_mass = pickle.load(f)
    K_BOL = pickle.load(f)
    core_tpd = pickle.load(f)


############### Plot Results
####### Define Set of Markers for Plotting
mrkrs = itertools.cycle(mrk.MarkerStyle.filled_markers)

###### Define Set of Line Styles for Plotting
lines = itertools.cycle(['-', '--', ':', '-.'])

####### BOL Keff vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = '3D Core BOL $K_{eff}$ vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.yaxis.set_major_locator(tck.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(tck.MultipleLocator(0.025))
ax.set_ylabel('BOL $K_{eff}$')
ax.grid()

for i in range(0, 7):
    # Plot BOL Keff Data with Error Bars for Each Moderator
    ax.errorbar(Rings, [j.n for j in K_BOL[i]], yerr=[j.s for j in K_BOL[i]],
                label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save 3D Reactor Core BOL Keff vs Fuel Element Radial Rings Plot
fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_BOL_Keff_Rings.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Reactor Core Radius vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure()
ax = fig.add_subplot(111)
str = '3D Core Radius vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_ylabel('3D Reactor Core Radius (cm)')
ax.grid()

for i in range(0, 7):
    # Plot Reactor Radius Data for Each Moderator
    ax.plot(Rings, core_rad[i],
            label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save 3D Reactor Core Radius vs Fuel Element Radial Rings Plot
ax.legend(loc='upper left')
fig.tight_layout()
fig.savefig(Pikname + '_Radius_Rings.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Reactor Core Height vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure()
ax = fig.add_subplot(111)
str = '3D Core Height vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_ylabel('3D Reactor Core Height (cm)')
ax.grid()

for i in range(0, 7):
    # Plot Reactor Height Data for Each Moderator
    ax.plot(Rings, core_hei[i],
            label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save 3D Reactor Core Height vs Fuel Element Radial Rings Plot
ax.legend(loc='upper left')
fig.tight_layout()
fig.savefig(Pikname + '_Height_Rings.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Reflector Mass vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure()
ax = fig.add_subplot(111)
str = 'Reflector Mass vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_yscale('log')
ax.set_ylabel('3D Reactor Radial Reflector Mass (Kg)')
ax.grid()

for i in range(0, 7):
    # Plot Radial Reflector Mass Data for Each Reflector
    ax.plot(Rings, ref_mass[i],
            label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save Reflector Mass vs Fuel Element Radial Rings Plot
ax.legend(loc='upper left')
fig.tight_layout()
fig.savefig(Pikname + '_RefMass_Rings.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Minimum Reactor Core Mass vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure()
ax = fig.add_subplot(111)
str = 'Minimum 3D Core Mass vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_yscale('log')
ax.set_ylabel('Minimum 3D Reactor Core Mass (Kg)')
ax.grid()

for i in range(0, 7):
    # Plot Reactor Mass Data for Each Moderator
    ax.plot(Rings, core_mass[i],
            label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save Minimum 3D Reactor Core Mass vs Fuel Element Radial Rings Plot
ax.legend(loc='upper left')
fig.tight_layout()
fig.savefig(Pikname + '_CoreMass_Rings.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Reactor Core Power Density vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = '3D Core Power Density vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_ylabel('3D Reactor Core Power Density (MW/m$^3$)')
ax.grid()

for i in range(0, 7):
    # Plot Reactor Power Density Data for Each Moderator
    ax.plot(Rings, core_tpd[i],
            label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save 3D Reactor Core Power Density vs Fuel Element Radial Rings Plot
fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_TPD_Rings.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()
