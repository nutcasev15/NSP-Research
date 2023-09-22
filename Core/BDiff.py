# OpenMC 0.13.3 Python Post Processing Script for Comparing Reflector BOL Keff Studies


############### Library Imports
import pickle
import itertools
import numpy as np
import matplotlib.markers as mrk
import matplotlib.pyplot as plt
import matplotlib.ticker as tck


############### Load Target Results from Pickle File
Tgtname = 'BOL-Nat-0'
with open(Tgtname + '.pckl', 'rb') as f:
    Rings_tgt = pickle.load(f)
    pow_tgt = pickle.load(f)
    Ref_thic_tgt = pickle.load(f)
    MF_Mod_tgt = pickle.load(f)
    Mname_tgt = pickle.load(f)
    Rname_tgt = pickle.load(f)
    core_rad_tgt = pickle.load(f)
    core_hei_tgt = pickle.load(f)
    ref_mass_tgt = pickle.load(f)
    core_mass_tgt = pickle.load(f)
    K_BOL_tgt = pickle.load(f)
    core_tpd_tgt = pickle.load(f)


############### Load Reference Results from Pickle File
Refname = 'BOL-Ref-0'
with open(Refname + '.pckl', 'rb') as f:
    Rings_ref = pickle.load(f)
    pow_ref = pickle.load(f)
    Ref_thic_ref = pickle.load(f)
    MF_Mod_ref = pickle.load(f)
    Mname_ref = pickle.load(f)
    Rname_ref = pickle.load(f)
    core_rad_ref = pickle.load(f)
    core_hei_ref = pickle.load(f)
    ref_mass_ref = pickle.load(f)
    core_mass_ref = pickle.load(f)
    K_BOL_ref = pickle.load(f)
    core_tpd_ref = pickle.load(f)


############### Find % Error From Reference Values
K_BOL_diff = np.divide((np.array(K_BOL_tgt) - np.array(K_BOL_ref)),
                       np.array(K_BOL_ref)) * 100
core_rad_diff = np.divide((np.array(core_rad_tgt) - np.array(core_rad_ref)),
                          np.array(core_rad_ref)) * 100
core_hei_diff = np.divide((np.array(core_hei_tgt) - np.array(core_hei_ref)),
                          np.array(core_hei_ref)) * 100
ref_mass_diff = np.divide((np.array(ref_mass_tgt) - np.array(ref_mass_ref)),
                           np.array(ref_mass_ref)) * 100
core_mass_diff = np.divide((np.array(core_mass_tgt) - np.array(core_mass_ref)),
                           np.array(core_mass_ref)) * 100
core_tpd_diff = np.divide((np.array(core_tpd_tgt) - np.array(core_tpd_ref)),
                          np.array(core_tpd_ref)) * 100


############### Plot Results
####### Define Set of Markers for Plotting
mrkrs = itertools.cycle(mrk.MarkerStyle.filled_markers)

###### Define Set of Line Styles for Plotting
lines = itertools.cycle(['-', '--', ':', '-.'])

####### BOL Keff % Difference vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = 'Difference in 3D Core BOL $K_{eff}$ vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_ylabel('Difference BOL $K_{eff}$ (%)')
ax.grid()

for i in range(0, 7):
    # Plot BOL Keff Difference Data with Error Bars for Each Moderator
    ax.errorbar(Rings_ref, [j.n for j in K_BOL_diff[i]],
                yerr=[j.s for j in K_BOL_diff[i]],
                label=Mname_ref[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save 3D Reactor Core BOL Keff % Difference Plot
fig.legend(loc='outside right upper')
fig.savefig(Tgtname + '_BOL_Keff_Rings_Diff.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Reactor Core Radius % Difference vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = 'Difference in 3D Core Radius vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_ylabel('Difference in 3D Reactor Core Radius (%)')
ax.grid()

for i in range(0, 7):
    # Plot Reactor Radius Difference Data for Each Moderator
    ax.plot(Rings_ref, core_rad_diff[i],
            label=Mname_ref[i][0], marker=next(mrkrs), linestyle=next(lines))


# Finalise and Save 3D Reactor Core Radius % Difference Plot
fig.legend(loc='outside right upper')
fig.savefig(Tgtname + '_Radius_Rings_Diff.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Reactor Core Height % Difference vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = 'Difference in 3D Core Height vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_ylabel('Difference in 3D Reactor Core Height (%)')
ax.grid()

for i in range(0, 7):
    # Plot Reactor Height Difference Data for Each Moderator
    ax.plot(Rings_ref, core_hei_diff[i],
            label=Mname_ref[i][0], marker=next(mrkrs), linestyle=next(lines))


# Finalise and Save 3D Reactor Core Height % Difference Plot
fig.legend(loc='outside right upper')
fig.savefig(Tgtname + '_Height_Rings_Diff.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Reflector Mass % Difference vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = 'Difference in Reflector Mass vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_ylabel('Difference in Radial Reflector Mass (%)')
ax.grid()

for i in range(0, 7):
    # Plot Reflector Mass Difference Data for Each Moderator
    ax.plot(Rings_ref, ref_mass_diff[i],
            label=Mname_ref[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save Reflector Mass % Difference Plot
fig.legend(loc='outside right upper')
fig.savefig(Tgtname + '_RefMass_Rings_Diff.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Minimum Reactor Core Mass % Difference vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = 'Difference in Minimum 3D Core Mass vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_ylabel('Difference in Minimum 3D Reactor Core Mass (%)')
ax.grid()

for i in range(0, 7):
    # Plot Reactor Mass Difference Data for Each Moderator
    ax.plot(Rings_ref, core_mass_diff[i],
            label=Mname_ref[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save Minimum 3D Reactor Core Mass % Difference Plot
fig.legend(loc='outside right upper')
fig.savefig(Tgtname + '_Mass_Rings_Diff.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Reactor Core Power Density % Difference vs Reactor Core Radial Rings for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
str = 'Difference in 3D Core Power Density vs Number of Radial Fuel Element Rings'
ax.set_title(str)
ax.set_xlabel('Fuel Element Radial Rings in Reactor Core')
ax.set_ylabel('Difference in 3D Reactor Core Thermal Power Density (%)')
ax.grid()

for i in range(0, 7):
    # Plot Reactor Power Density Difference Data for Each Moderator
    ax.plot(Rings_ref, core_tpd_diff[i],
            label=Mname_ref[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save 3D Reactor Core Power Density % Difference Plot
fig.legend(loc='outside right upper')
fig.savefig(Tgtname + '_TPD_Rings_Diff.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()
