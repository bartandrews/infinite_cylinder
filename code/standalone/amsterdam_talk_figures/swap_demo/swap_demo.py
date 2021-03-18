import numpy as np
import matplotlib
import csv
import matplotlib.pyplot as plt
import sys
from itertools import product
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
import matplotlib.patches as patches
from numpy.polynomial.polynomial import polyfit
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import ConnectionPatch
from scipy import stats
from fractions import Fraction


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    fig = plt.figure(figsize=(10, 1))

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    ax1 = plt.subplot(111)

    ax1.plot(0, 0, 'ko', markerfacecolor='white')
    ax1.plot(1, 0, 'ko')
    ax1.plot(2, 0, 'ko')
    ax1.plot(3, 0, 'ko', markerfacecolor='white')
    ax1.plot(4, 0, 'ko', markerfacecolor='white')
    ax1.plot(5, 0, 'ko')
    ax1.plot(6, 0, 'ko')
    ax1.plot(7, 0, 'ko', markerfacecolor='white')

    ax1.annotate("0", xy=(0, 0), xytext=(0, 5), textcoords='offset points', fontsize=12)
    ax1.annotate("1", xy=(1, 0), xytext=(0, 5), textcoords='offset points', fontsize=12)
    ax1.annotate("5", xy=(2, 0), xytext=(0, 5), textcoords='offset points', fontsize=12)
    ax1.annotate("4", xy=(3, 0), xytext=(0, 5), textcoords='offset points', fontsize=12)
    ax1.annotate("3", xy=(4, 0), xytext=(0, 5), textcoords='offset points', fontsize=12)
    ax1.annotate("2", xy=(5, 0), xytext=(0, 5), textcoords='offset points', fontsize=12)
    ax1.annotate("6", xy=(6, 0), xytext=(0, 5), textcoords='offset points', fontsize=12)
    ax1.annotate("7", xy=(7, 0), xytext=(0, 5), textcoords='offset points', fontsize=12)

    ax1.axvspan(1.5, 5.5, alpha=0.25, color='green')
    ax1.hlines(y=0, xmin=0, xmax=8, linewidth=2, color='k')
    ax1.hlines(y=0, xmin=-1, xmax=0, linewidth=2, color='k', linestyles='dashed')

    plt.axis('off')

    ####################################################################################################################

    plt.savefig("/home/bart/Documents/presentations/2021_03_18/figures/swap_demo_6.png", bbox_inches='tight', dpi=300)
    plt.show()
