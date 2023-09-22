# OpenMC 0.13.3 Python Post Processing Script for MF Ratio Depletion Study


############### Library Imports
import pickle
import itertools
import numpy as np
import matplotlib.markers as mrk
import matplotlib.pyplot as plt
import matplotlib.ticker as tck


############### Load Results from Pickle File
Pikname = 'Dep-Ref'
with open(Pikname + '.pckl', 'rb') as f:
    MF_vol = pickle.load(f)
    dep_sch = pickle.load(f)
    Mname = pickle.load(f)
    MF_den_ratio = pickle.load(f)
    chnl_dia = pickle.load(f)
    hex_dia = pickle.load(f)
    hex_mass = pickle.load(f)
    K_EVL = pickle.load(f)
    K_BOL = pickle.load(f)
    K_EOL = pickle.load(f)
    U_burn = pickle.load(f)


############### Plot Results
####### Define Set of Markers for Plotting
mrkrs = itertools.cycle(mrk.MarkerStyle.filled_markers)

###### Define Set of Line Styles for Plotting
lines = itertools.cycle(['-', '--', ':', '-.'])

####### BOL Kinf vs MF_vol for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
ax.set_title('BOL $K_{inf}$ vs Moderator to Fuel Volume Ratio')
ax.set_xscale('log')
ax.set_xlabel('Moderator to Fuel Volume Ratio')
ax.autoscale(False, axis='y')
ax.set_ylim(0.6, 1.8)
ax.yaxis.set_major_locator(tck.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(tck.MultipleLocator(0.025))
ax.set_ylabel('BOL $K_{inf}$')
ax.grid()

for i in range(0, 7):
    # Plot BOL Kinf Data with Error Bars for Each Moderator
    ax.errorbar(MF_vol, [j[0] for j in K_BOL[i]], yerr=[j[1] for j in K_BOL[i]],
                label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save BOL Kinf vs MF_vol Plot
fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_BOL_MFvol.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()


####### EOL Kinf vs MF_vol for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
ax.set_title('EOL $K_{inf}$ vs Moderator to Fuel Volume Ratio')
ax.set_xscale('log')
ax.set_xlabel('Moderator to Fuel Volume Ratio')
ax.autoscale(False, axis='y')
ax.set_ylim(0.6, 1.8)
ax.yaxis.set_major_locator(tck.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(tck.MultipleLocator(0.025))
ax.set_ylabel('EOL $K_{inf}$')
ax.grid()

for i in range(0, 7):
    # Plot EOL Kinf Data with Error Bars for Each Moderator
    ax.errorbar(MF_vol, [j[0] for j in K_EOL[i]], yerr=[j[1] for j in K_EOL[i]],
                label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))


# Finalise and Save EOL Kinf vs MF_vol Plot
fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_EOL_MFvol.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Difference in Kinf vs MF_vol for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
ax.set_title('Difference in $K_{inf}$ vs Moderator to Fuel Volume Ratio')
ax.set_xscale('log')
ax.set_xlabel('Moderator to Fuel Volume Ratio')
ax.autoscale(False, axis='y')
ax.set_ylim(-0.3, 0)
ax.yaxis.set_major_locator(tck.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(tck.MultipleLocator(0.025))
ax.set_ylabel('Difference in $K_{inf}$')
ax.grid()

for i in range(0, 7):
    # Calculate Difference in Kinf between BOL and EOL
    K_diff = (np.array(K_EOL) - np.array(K_BOL)).tolist()
    # Plot Kinf Difference Data for Each Moderator
    ax.plot(MF_vol, [j[0] for j in K_diff[i]],
            label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))


# Finalise and Save Kinf Difference vs MF_vol Plot
fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_KDiff_MFvol.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### U Atom Burnup at EOL vs MF_vol for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
ax.set_title('U Atom Burnup vs Moderator to Fuel Volume Ratio')
ax.set_xscale('log')
ax.set_xlabel('Moderator to Fuel Volume Ratio')
ax.autoscale(False, axis='y')
ax.set_ylim(4, 8)
ax.yaxis.set_major_locator(tck.MultipleLocator(1))
ax.yaxis.set_minor_locator(tck.MultipleLocator(0.25))
ax.set_ylabel('U Burnup at EOL (Atom %)')
ax.grid()

for i in range(0, 7):
    # Plot U Burnup Data for Each Moderator
    ax.plot(MF_vol, U_burn[i], label=Mname[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save U Atom % Burnup vs MF_vol
fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_Burn_MFvol.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()
