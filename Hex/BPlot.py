# OpenMC 0.15.0 Python Post Processing Script for BOL Kinf Study


############### Library Imports
import pickle
import itertools
import numpy as np
import matplotlib.markers as mrk
import matplotlib.pyplot as plt
import matplotlib.ticker as tck


############### Load Results from Pickle File
PklName = 'BOL-Ref'
with open(PklName + '.pkl', 'rb') as f:
    MF_vol = pickle.load(f)
    ModName = pickle.load(f)
    MF_den_ratio = pickle.load(f)
    channel_dia = pickle.load(f)
    hex_dia = pickle.load(f)
    hex_mass = pickle.load(f)
    K_BOL = pickle.load(f)


############### Plot Results
####### Define Set of Markers for Plotting
marks = itertools.cycle(mrk.MarkerStyle.filled_markers)

###### Define Set of Line Styles for Plotting
lines = itertools.cycle(['-', '--', ':', '-.'])

####### BOL Kinf vs MF_vol for Various Moderators
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('BOL $K_{inf}$ vs Moderator to Fuel Volume Ratio')
ax.set_xlim(left=1, right=10)
ax.set_xlabel('Moderator to Fuel Volume Ratio')
ax.yaxis.set_major_locator(tck.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(tck.MultipleLocator(0.025))
ax.set_ylabel('BOL $K_{inf}$')
ax.grid()

for Kinf, MName in zip(K_BOL, ModName):
    # Plot BOL Kinf Data with Error Bars for Each Moderator
    ax.errorbar(MF_vol, [j.n for j in Kinf], yerr=[j.s for j in Kinf],
                label=MName[0], marker=next(marks), linestyle=next(lines))

# Finalise and Save BOL Kinf vs MF_vol Plot
ax.legend(loc='lower left')
fig.tight_layout()
fig.savefig(PklName + '_BOL_Kinf_MFvol.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()


####### Optimum MF_vol for Various Moderators
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Optimum Moderator to Fuel Volume Ratios')
ax.tick_params('x', rotation=30, labelsize=10)
ax.set_ylabel('Optimum MF Volume Ratio')

# Define Empty Container to Store Optimum MF_vol
MF_optima = []

for Kinf in K_BOL:
    # Fit BOL Kinf Data to 4th Order Polynomial
    p = np.polynomial.Polynomial.fit(MF_vol, [j.n for j in Kinf], 4)

    # Find Maximum by Evaluating the Polynomial for a Large Number of Points
    MF = np.linspace(1, 10, 512)
    K = p(MF)
    MF_optima.append(MF[K.argmax()])

# Plot Optimum MF_vol Data for Each Moderator
ax.bar_label(ax.bar([i[0] for i in ModName], MF_optima, width=0.75),
             fmt='%.3f')

# Finalise and Save Optimum MF_vol Plot
fig.tight_layout()
fig.savefig(PklName + '_Optima_MFvol.pdf', format='pdf')

# Clear Figure, Axes and Container
fig.clear()
ax.clear()


####### Optimum MF_mass for Various Moderators
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Optimum Moderator to Fuel Mass Ratios')
ax.tick_params('x', rotation=30, labelsize=10)
ax.set_ylabel('Optimum MF Mass Ratio')

# Define Empty Container to Store Optimum MF_mass
MF_optima = []

for Kinf, ModDens in zip(K_BOL, MF_den_ratio):
    # Fit BOL Kinf Data to 4th Order Polynomial
    p = np.polynomial.Polynomial.fit(MF_vol, [j.n for j in Kinf], 4)

    # Find Maximum by Evaluating the Polynomial for a Large Number of Points
    MF = np.linspace(1, 10, 512)
    K = p(MF)
    MF_optima.append(float(MF[K.argmax()]) * ModDens[0])

# Plot Optimum MF_mass Data for Each Moderator
ax.bar_label(ax.bar([i[0] for i in ModName], MF_optima, width=0.75),
             fmt='%.3f')

# Finalise and Save Optimum MF_mass Plot
fig.tight_layout()
fig.savefig(PklName + '_Optima_MFmass.pdf', format='pdf')

# Clear Figure, Axes and Container
fig.clear()
ax.clear()


####### Fuel Element Linear Mass vs MF_vol for Various Moderators
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('FE Linear Mass vs Moderator to Fuel Volume Ratio')
ax.set_xlim(left=1, right=10)
ax.set_xlabel('Moderator to Fuel Volume Ratio')
ax.set_ylabel('Fuel Element Linear Mass (kg/m)')
ax.grid()

for i in range(0, 8):
    # Plot Linear Mass Data for Each Moderator
    ax.plot(MF_vol, hex_mass[i],
            label=ModName[i][0], marker=next(marks), linestyle=next(lines))

# Finalise and Save FE Linear Mass vs MF_vol Plot
ax.legend(loc='upper left')
fig.tight_layout()
fig.savefig(PklName + '_Mass_MFvol.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()
