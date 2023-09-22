# OpenMC 0.13.3 Python Post Processing Script for 3D Core Reflector Depletion Study


############### Library Imports
import pickle
import itertools
import matplotlib.markers as mrk
import matplotlib.pyplot as plt
import matplotlib.ticker as tck


############### Load Results from Pickle File
Pikname = 'Dep-Ref-0'
with open(Pikname + '.pckl', 'rb') as f:
    Rings = pickle.load(f)
    pow = pickle.load(f)
    Ref_thic = pickle.load(f)
    dep_sch = pickle.load(f)
    MF_Mod = pickle.load(f)
    Mname = pickle.load(f)
    Rname = pickle.load(f)
    core_rad = pickle.load(f)
    core_hei = pickle.load(f)
    ref_mass = pickle.load(f)
    core_mass = pickle.load(f)
    K_EVL = pickle.load(f)
    K_BOL = pickle.load(f)
    K_EOL = pickle.load(f)
    U_burn = pickle.load(f)


############### Plot Results
####### Define Set of Markers for Plotting
mrkrs = itertools.cycle(mrk.MarkerStyle.filled_markers)

###### Define Set of Line Styles for Plotting
lines = itertools.cycle(['-', '--', ':', '-.'])

####### 3D Reactor Core BOL Keff vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = '3D Core BOL $K_{eff}$ vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.autoscale(False, axis='y')
ax.set_ylim(0.8, 1.8)
ax.yaxis.set_major_locator(tck.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(tck.MultipleLocator(0.025))
ax.set_ylabel('BOL $K_{eff}$')
ax.grid()

for i in range(0, 5):
    # Plot 3D Core BOL Keff Data with Error Bars for Each Moderator
    ax.errorbar(Rings, [j[0] for j in K_BOL[i]], yerr=[j[1] for j in K_BOL[i]],
                label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save 3D Core BOL Keff vs Fuel Element Radial Rings Plot
fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_Rings_BOL.png', format='png')

# Clear Figure and Axes
fig.clear()
ax.clear()


####### EOL Keff vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = '3D Core EOL $K_{eff}$ vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.autoscale(False, axis='y')
ax.set_ylim(0.8, 1.8)
ax.yaxis.set_major_locator(tck.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(tck.MultipleLocator(0.025))
ax.set_ylabel('EOL $K_{eff}$')
ax.grid()

for i in range(0, 5):
    # Plot 3D Core EOL Keff Data with Error Bars for Each Moderator
    ax.errorbar(Rings, [j[0] for j in K_EOL[i]], yerr=[j[1] for j in K_EOL[i]],
                label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save 3D Core EOL Keff vs Fuel Element Radial Rings Plot
fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_Rings_EOL.png', format='png')

# Clear Figure and Axes
fig.clear()
ax.clear()


####### U Atom Burnup at EOL vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = 'U Atom Burnup vs Number of Radials Fuel Elements Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.autoscale(False, axis='y')
ax.set_ylim(0, 8)
ax.yaxis.set_major_locator(tck.MultipleLocator(1))
ax.yaxis.set_minor_locator(tck.MultipleLocator(0.25))
ax.set_ylabel('U Burnup at EOL (Atom %)')
ax.grid()

for i in range(0, 5):
    # Plot U Burnup Data for Each Moderator
    ax.plot(Rings, U_burn[i], label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save U Atom % Burnup vs Fuel Element Radial Rings Plot
fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_Rings_Burn.png', format='png')

# Clear Figure and Axes
fig.clear()
ax.clear()
