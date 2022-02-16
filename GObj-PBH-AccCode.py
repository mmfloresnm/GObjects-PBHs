#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GObj-PBH-AccCode.py -- Generates accretion plot seen in 
                       G Objects and Primordial Black Holes by M. Flores, 
                       A. Kusenko, A. M. Ghez & S. Naoz

For complete details see article.

Written by Marcos M. Flores
Last Update: 11, Feb. 2022
"""


import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


def MEdd(Mm):
    """
    MEdd(Mm) - Defines the Eddington mass rate as where Mm = M/M_\odot, i.e.
               the mass normalized by a solar-mass
    """
    return 1.402e17 * Mm

def MBondi(Mm, Nn, Tt):
    """
    MBondi(Mm, Nn, Tt) - Defines the Bondi accretion mass rate as specified by 
                         eq. 11. Here, Mm is the same as in the above 
                         funciton, Nn is the normalized number densty with 
                         Nn = n/(1 cm^-3) and Tt is the normalized 
                         temperature for Tt = T/(10^4 K)
    """
    return 2.4e10 * Mm**2 * Nn * Tt**(-3/2)

def mdot(Mm, Nn, Tt):
    """
    mdot(Mm, Nn, Tt) - Defines the dimensionless acrretion rate
    """
    return MBondi(Mm, Nn, Tt)/ MEdd(Mm)

npts = int(1e3) # Number of points for plotting

mdots = np.logspace(-7, 4, npts) # Define mdot axis range

lumeta = lambda eta: eta*mdots**2 # Dimensionless luminosity as function of
                                  # \eta and \mdot
lumeps = lambda eps: eps*mdots    # Dimensionless luminosity as a function of
                                  # effeciency \epsilon and \mdot

# Plotting Information

fig, ax = plt.subplots()

# Data from Park & Ostriker, ApJ 549, 100 (2001)
WarmSphData = np.array([[5.068942952780017, 0.0004993091965674814], 
[9.971583362138187, 0.000641305873355304],
[20.13770948421302, 0.0012370615382025437],
[30.14881841137572, 0.0020090425294644926],
[39.63014573197902, 0.003162277660168376], 
[51.41309556502552, 0.0039364858588476544],
[99.9240254677894, 0.025723264710910673],
[80.06600696100007, 0.012526394555458566],
[71.21098263848482,0.01006276640779734], 
[60.11536761703846, 0.00629380152731093]])

ColdSphData = np.array([[178.68493523769783, 0.000015495920916907527], 
[70.92834446803616, 8.823729256112076e-6], 
[42.698414991159154, 6.658395148989865e-6], 
[28.522652468547943,  4.794093741105753e-6],
[7.178535174035825, 1.5065693591480495e-6],
[4.265089513340143, 9.276652840859881e-7],
[2.8484111522566127, 4.378288223191186e-7],
[0.7073866701120404, 7.713087547747421e-8], 
[0.3594840740082596, 3.528209360426702e-8], 
[0.1390041176046312, 1.4240741682392281e-8], 
[0.11889560658639713, 1.0915496086248792e-8]])

HotSphData = np.array([[0.0010137567026216356, 3.142552983956818e-11], 
[0.010031334112568163, 2.5087450104207354e-9], 
[0.01974123388158432, 6.4130587335530395e-9],
[0.039868318305332794, 1.27637491673712e-8],
[0.060460622055860565, 1.5888648468825485e-8],
[0.07838566985889764, 6.215531343872624e-9], 
[0.0916290711690162, 6.215531343872624e-9]])

# Data plots

WarmPlt, = ax.loglog(WarmSphData[:,0], WarmSphData[:,1], "o", color="purple",
            label=r"${\rm Warm\ Spherical}$")

ColdPlt, = ax.loglog(ColdSphData[:,0], ColdSphData[:,1], "s", color="darkblue",
            label= r"${\rm Cold\ Spherical}$")

HotPlt, = ax.loglog(HotSphData[:,0], HotSphData[:,1], "^", color="darkred",
            label=r"${\rm Hot\ Spherical}$")

ax.legend(handles=[WarmPlt, ColdPlt, HotPlt], prop={'size': 13})

# Flat lines
ax.loglog(mdots, mdots*0 + 1*1e-3, "b--")
ax.loglog(mdots, mdots*0 + 1, color="gray", linestyle="--")

ax.fill_between(mdots, lumeta(1.0e-5), lumeta(1), color='lightsteelblue',
                alpha=0.4)

# Lines of const. \epsilon
ax.loglog(mdots, lumeps(1), "k--")
ax.loglog(mdots, lumeps(1.0e-5), "k--")

# Test mdot values

TtEx = 550/(1e4)
NnMin = 0.2e5
NnMax = 2e5

ax.axvspan(mdot(1, NnMin, TtEx), mdot(1, NnMax, TtEx),color='green',
           alpha=0.25)


# Luminosity lables
ax.text(1.25e-6, 1.5e-3, r"$\mathrm{Inferred\ Luminosity}$", fontsize=14)
ax.text(1.25e-6, 1.5e0, r"$\mathrm{Eddington\ Luminosity}$", fontsize=14)

# Epsilon labels
ax.text(2e-6, 5e-7, r"$\epsilon = 1$", rotation=30, fontsize=16)
ax.text(2e-6, 5e-11, r"$\epsilon = 10^{-5}$", rotation=30, fontsize=16)

# Axis adjustments

ax.set_xlabel(r"$\dot{m}$", fontsize=20)
ax.set_ylabel(r"$l$", fontsize=20)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

ax.set_xlim(1e-6, 1.0e4)
ax.set_ylim(1e-12, 10)

# Set figure size
fig.set_size_inches(8.0,6.0)