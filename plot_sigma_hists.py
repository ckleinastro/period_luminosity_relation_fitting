from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from numpy import array, log10, loadtxt, exp
from scipy import interpolate

import sys

plot_dir = "plots/"

band_list =              ["U",  "B", "hipp",  "V",  "R",  "I",  "z",  "J", "H", "K", "W1", "W2", "W3"]

data_dir = "/Volumes/Extra_HDD/Research/Multiband_PLR_Fitting/lowexpandedEBV_trace_data/"
plr_db_prefix = "M_lowexpandedEBV"
# Removed "7" and "9" because they don't converge until well into the analysis period of the trace
trace_run_nums = [1, 2, 3, 4, 5, 6, 8] 

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)


# for band in band_list:
for band in ["hipp"]:

    full_sigma_trace = []
    for n in trace_run_nums:
        sigma_trace = loadtxt(data_dir + band + "_" + plr_db_prefix + "_1" + str(n) + "_sigma_trace.txt")
        full_sigma_trace.append(sigma_trace.tolist())
        
    full_sigma_trace = array(full_sigma_trace)
    full_sigma_trace = full_sigma_trace.flatten()



    fig = plt.figure(figsize=(3.3, 3.3))
    ax1 = subplot(111)
    
    ax1.hist(full_sigma_trace, bins=100, histtype="stepfilled", color="gray")

    ax1.set_xlabel(r"Intrinsic Scatter")
    ax1.set_ylabel(r"Count")

    ax1.set_xlim(0, full_sigma_trace.max())

    # This code draws major and minor tick lines. Major ticks get number labels.
    majorLocator_x = MultipleLocator(0.02)
    minorLocator_x = MultipleLocator(0.01)
    ax1.xaxis.set_major_locator(majorLocator_x)
    ax1.xaxis.set_minor_locator(minorLocator_x)

    majorLocator_y1 = MultipleLocator(2500)
    minorLocator_y1 = MultipleLocator(1250)
    ax1.yaxis.set_major_locator(majorLocator_y1)
    ax1.yaxis.set_minor_locator(minorLocator_y1)

    ax1.text(0.88, 0.92, r"$%s$" % band, transform=ax1.transAxes, fontsize=10)

    # pos =         [left, bottom, width, height]
    ax1.set_position([0.22, 0.14, 0.75, 0.82])

    canvas = FigureCanvas(fig)
    canvas.print_figure(plot_dir + band + "_sigma.pdf", dpi=300)
    close("all")
