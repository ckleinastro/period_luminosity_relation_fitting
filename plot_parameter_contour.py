from pylab import *
from matplotlib.pylab import *
from matplotlib.ticker import NullFormatter
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from numpy import array, log10, loadtxt
from scipy import interpolate

import sys

plot_dir = "plots/"

band_list =              ["U",  "B", "hipp",  "V",  "R",  "I",  "z",  "J", "H", "K", "W1", "W2", "W3"]


P_0 = 0.52853966619770265

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






from numpy import *
from scipy.interpolate import interp2d
import scipy.stats as stats

def percent_level(z, low_val):
    n_pix_total = float(z.shape[0]*z.shape[1])
    n_pix_above = float(where(z>low_val)[0].shape[0])
    return float(n_pix_above/n_pix_total)

def find_level(z, target_percent):
    n_pix_total = float(z.shape[0]*z.shape[1])
    new_level = z.mean()/(n_pix_total)
    current_percent_level = percent_level(z, new_level)
    while abs(current_percent_level - target_percent) > 0.001:
        if current_percent_level < target_percent:
            new_level = new_level*0.99
        else:
            new_level = new_level*1.01
        current_percent_level = percent_level(z, new_level)
    return new_level






# for band in band_list:
for band in ["W3"]:
    full_M_0_trace = []
    full_alpha_trace = []
    for n in trace_run_nums:
    # for n in [1]:
        M_0_trace = loadtxt(data_dir + band + "_" + plr_db_prefix + "_1" + str(n) + "_M_0_trace.txt")
        alpha_trace = loadtxt(data_dir + band + "_" + plr_db_prefix + "_1" + str(n) + "_alpha_trace.txt")
        full_M_0_trace.append(M_0_trace.tolist())
        full_alpha_trace.append(alpha_trace.tolist())
        
    full_M_0_trace = array(full_M_0_trace)
    full_M_0_trace = full_M_0_trace.flatten()
    full_alpha_trace = array(full_alpha_trace)
    full_alpha_trace = full_alpha_trace.flatten()

    x_min = full_M_0_trace.mean()-5*full_M_0_trace.std()
    x_max = full_M_0_trace.mean()+5*full_M_0_trace.std()
    y_min = full_alpha_trace.mean()-5*full_alpha_trace.std()
    y_max = full_alpha_trace.mean()+5*full_alpha_trace.std()
    
    x_data = array([full_M_0_trace]).T
    y_data = array([full_alpha_trace]).T
    rvs = np.append(x_data, y_data, axis=1)

    kde = stats.kde.gaussian_kde(rvs.T)

    # Regular grid to evaluate kde upon
    x_flat = np.r_[x_min:x_max:256j]
    y_flat = np.r_[y_min:y_max:256j]
    x,y = np.meshgrid(x_flat,y_flat)
    grid_coords = np.append(x.reshape(-1,1),y.reshape(-1,1),axis=1)

    z = kde(grid_coords.T)
    z = z.reshape(256,256)



    l_999 = find_level(z, 0.001)
    l_99 = find_level(z, 0.01)
    l_975 = find_level(z, 0.025)
    l_95 = find_level(z, 0.05)
    l_90 = find_level(z, 0.10)
    l_85 = find_level(z, 0.15)
    l_80 = find_level(z, 0.20)
    l_70 = find_level(z, 0.30)
    levels_list=[l_70, l_80, l_85, l_90, l_95, l_975, l_99, l_999, z.max()]
    
    
    
    nullfmt = NullFormatter()   # no labels
    
    # definitions for the axes
    left, width = 0.17, 0.58
    bottom, height = 0.17, 0.58
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    
    rect_band_label = [left_h, bottom_h, 0.2, 0.2]
    
    fig = plt.figure(figsize=(3.3, 3.3))
    
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    axBandLabel = plt.axes(rect_band_label)
    
    axBandLabel.yaxis.set_ticks([])
    axBandLabel.xaxis.set_ticks([])
    axBandLabel.set_xlim(0,1)
    axBandLabel.set_ylim(0,1)
    axBandLabel.text(0.5, 0.5, r"$%s$" % band, fontsize=18, ha='center', va='center')
    

    
    axScatter.contourf(x, y, z, levels=levels_list, origin="lower", cmap=cm.gist_yarg)
    axScatter.contour(x, y, z, levels=levels_list, origin="lower", colors='k')
    # ax1.scatter(x_data, y_data, alpha=0.15, color="blue", linewidths=0)
    axScatter.errorbar(x_data.mean(), y_data.mean(), y_data.std(), x_data.std(), marker="o", color="r")
    axScatter.set_xlabel(r"$M_0$")
    axScatter.set_ylabel(r"$\alpha$")
    axScatter.set_xlim(x_min, x_max)
    axScatter.set_ylim(y_min, y_max)
    
    
    # This code draws major and minor tick lines. Major ticks get number labels.
    majorLocator_x = MultipleLocator(0.03)
    minorLocator_x = MultipleLocator(0.015)
    axScatter.xaxis.set_major_locator(majorLocator_x)
    axScatter.xaxis.set_minor_locator(minorLocator_x)

    majorLocator_y1 = MultipleLocator(0.3)
    minorLocator_y1 = MultipleLocator(0.15)
    axScatter.yaxis.set_major_locator(majorLocator_y1)
    axScatter.yaxis.set_minor_locator(minorLocator_y1)
    
    
    
    
    axHistx.hist(full_M_0_trace, bins=100, normed=False, histtype="stepfilled", color="gray", alpha=1.0)
    
    axHisty.hist(full_alpha_trace, bins=100, normed=False, histtype="stepfilled", orientation="horizontal", color="gray", alpha=1.0)
    
    axHistx.set_ylabel("Count")
    axHisty.set_xlabel("Count")
    
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )
    
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # This code draws major and minor tick lines. Major ticks get number labels.
    # hist_majorLocator_y = MultipleLocator(20)
    # hist_minorLocator_y = MultipleLocator(10)
    # axHistx.yaxis.set_major_locator(hist_majorLocator_y)
    # axHistx.yaxis.set_minor_locator(hist_minorLocator_y)

    # hist_majorLocator_x1 = MultipleLocator(1.0)
    # hist_minorLocator_x1 = MultipleLocator(0.5)
    # axHisty.xaxis.set_major_locator(hist_majorLocator_x1)
    # axHisty.xaxis.set_minor_locator(hist_minorLocator_x1)
    

    
    xticksHisty = axHisty.xaxis.get_major_ticks()
    for x_tick_label in xticksHisty:
        x_tick_label.label1.set_visible(False)
    
    yticksHistx = axHistx.yaxis.get_major_ticks()
    for y_tick_label in yticksHistx:
        y_tick_label.label1.set_visible(False)


    canvas = FigureCanvas(fig)
    canvas.print_figure(plot_dir + band + "_contour.pdf", dpi=300)
    close("all")




