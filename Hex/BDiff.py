# OpenMC 0.13.3 Python Post Processing Script for Comparing BOL Kinf Studies


############### Library Imports
import pickle
import itertools
import numpy as np
import matplotlib.markers as mrk
import matplotlib.pyplot as plt


############### Load Target Results from Pickle File
Tgtname = 'BOL-Nat'
with open(Tgtname + '.pckl', 'rb') as f:
    MF_vol_tgt = pickle.load(f)
    Mname_tgt = pickle.load(f)
    MF_den_ratio_tgt = pickle.load(f)
    chnl_dia_tgt = pickle.load(f)
    hex_dia_tgt = pickle.load(f)
    hex_mass_tgt = pickle.load(f)
    K_BOL_tgt = pickle.load(f)


############### Load Reference Results from Pickle File
Refname = 'BOL-Ref'
with open(Refname + '.pckl', 'rb') as f:
    MF_vol_ref = pickle.load(f)
    Mname_ref = pickle.load(f)
    MF_den_ratio_ref = pickle.load(f)
    chnl_dia_ref = pickle.load(f)
    hex_dia_ref = pickle.load(f)
    hex_mass_ref = pickle.load(f)
    K_BOL_ref = pickle.load(f)


############### Find % Error From Reference Values
K_BOL_diff = np.divide((np.array(K_BOL_tgt) - np.array(K_BOL_ref)),
                       np.array(K_BOL_ref)) * 100


############### Plot Results
####### Define Set of Markers for Plotting
mrkrs = itertools.cycle(mrk.MarkerStyle.filled_markers)

###### Define Set of Line Styles for Plotting
lines = itertools.cycle(['-', '--', ':', '-.'])

####### BOL Kinf % Difference vs MF_vol for Various Moderators
fig = plt.figure(layout='constrained')
ax = fig.add_subplot(111)
ax.set_title('Difference in BOL $K_{inf}$ vs Moderator to Fuel Volume Ratio')
ax.set_xscale('log')
ax.set_xlabel('Moderator to Fuel Volume Ratio')
ax.set_ylabel('Difference in BOL $K_{inf}$ (%)')
ax.grid()

for i in range(0, 7):
    # Plot BOL Kinf Difference Data with Relative Error for Each Moderator
    ax.errorbar(MF_vol_ref, [j.n for j in K_BOL_diff[i]], yerr=[j.s for j in K_BOL_diff[i]],
                label=Mname_ref[i][0], marker=next(mrkrs), linestyle=next(lines))

# Finalise and Save BOL Kinf % Difference vs MF_vol Plot
fig.legend(loc='outside right upper')
fig.savefig(Tgtname + '_BOL_Kinf_Diff.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

####### Optimum MF % Difference for Various Moderators
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Difference in Optimum Moderator to Fuel Ratios')
ax.tick_params('x', rotation=30, labelsize=10)
ax.set_ylabel('Difference in Optimum MF Ratio (%)')

# Define Empty Container to Store Optimum MF Difference
MF_optim_diff = []

for i in range(0, 7):
    # Fit BOL Kinf Data to 4th Order Polynomial
    p_ref = np.polynomial.Polynomial.fit(MF_vol_ref, [j.n for j in K_BOL_ref[i]], 4)
    p_tgt = np.polynomial.Polynomial.fit(MF_vol_tgt, [j.n for j in K_BOL_tgt[i]], 4)

    # Find Maximum by Evaluating the Polynomial for a Large Number of Points
    MF = np.logspace(0, np.log10(25), 512)
    K_tgt = p_tgt(MF)
    K_ref = p_ref(MF)
    MF_optim_diff.append(((MF[K_tgt.argmax()] - MF[K_ref.argmax()])
                          / MF[K_ref.argmax()]) * 100)

# Plot Optimum MF Ratio Difference Data for Each Moderator
ax.bar_label(ax.bar([i[0] for i in Mname_ref], MF_optim_diff, width=0.75),
             fmt='%.3f')

# Finalise and Save Optimum MF % Difference Ratio Plot
fig.tight_layout()
fig.savefig(Tgtname + '_Optim_Diff.pdf', format='pdf')

# Clear Figure, Axes and Container
fig.clear()
ax.clear()
del MF_optim_diff
