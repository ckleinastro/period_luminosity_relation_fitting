from numpy import loadtxt, linspace, array
from scipy import interpolate
from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib

# plr_db_prefix = "M_simpleEBV"
# trace_run_nums = [1, 2, 3, 4, 5, 6, 7, 9] # Remove "8" because it doesn't converge until well into the analysis period of the trace

plr_db_prefix = "M_lowexpandedEBV"
trace_run_nums = [1, 2, 3, 4, 5, 6, 8] # Remove "7" and "9" because they don't converge until well into the analysis period of the trace
n_trace_runs = 9
n_samples_analysis = 50000

data_dir = "/Volumes/Extra_HDD/Research/Multiband_PLR_Fitting/trace_data/"
band_list = ["U", "B", "hipp", "V", "R", "I", "z", "J", "H", "K", "W1", "W2", "W3"]



rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rc('axes', labelsize=10)






raw_red_color = loadtxt("red_mapping.txt")
raw_green_color = loadtxt("green_mapping.txt")
raw_blue_color = loadtxt("blue_mapping.txt")

red_function = interpolate.interp1d(raw_red_color[:,0], raw_red_color[:,1], kind="linear")
green_function = interpolate.interp1d(raw_green_color[:,0], raw_green_color[:,1], kind="linear")
blue_function = interpolate.interp1d(raw_blue_color[:,0], raw_blue_color[:,1], kind="linear")

def make_color(val):
    val *= 255
    color = [red_function(val)/255., green_function(val)/255., blue_function(val)/255., 1.0]
    return color

color_vals = linspace(1, 0, n_trace_runs)
colors_list = []
for cv in color_vals:
    colors_list.append(make_color(cv))





fig = plt.figure(figsize=(8, 12))
for m in range(len(band_list)):
    ax = fig.add_subplot(len(band_list),1,1+m)
    big_sigma_trace = []
    for n in trace_run_nums:
        sigma_trace = loadtxt(data_dir + band_list[m] + "_" + plr_db_prefix + "_1" + str(n) + "_sigma_trace.txt")[-n_samples_analysis:]
                
        big_sigma_trace.append(sigma_trace.tolist())
        ax.plot(sigma_trace, color=colors_list[n-1], alpha=0.1)
    big_sigma_trace = array(big_sigma_trace)
    big_sigma_mean = big_sigma_trace.mean()
    big_sigma_std = big_sigma_trace.std()
    ax.hlines(big_sigma_mean, 0, n_samples_analysis, colors="k")
    ax.hlines([big_sigma_mean+big_sigma_std, big_sigma_mean-big_sigma_std], 0, n_samples_analysis, colors="k", linestyles="--")
    
    ylabel_obj = ax.set_ylabel(band_list[m] + "\n" + r" %.4f$\pm$%.4f" % (big_sigma_mean, big_sigma_std))
    ylabel_obj.set_rotation(0)
    ylabel_obj.set_horizontalalignment("right")
    ax.set_xlim(0, n_samples_analysis)
    # pos = [left, bottom, width, height]
    ax_pos = ax.get_position().bounds
    ax.set_position([ax_pos[0]+0.04, ax_pos[1]-0.0105*(m-3.8), ax_pos[2]+0.03, ax_pos[3]+0.015])
    
    if m < len(band_list)-1:
        for tl in ax.get_xticklabels():
            tl.set_visible(False)
fig.suptitle("sigma traces for each band")
canvas = FigureCanvas(fig)
canvas.print_figure(plr_db_prefix + "_sigma_traces.png", dpi=300)
close("all")


fig = plt.figure(figsize=(8, 12))
for m in range(len(band_list)):
    ax = fig.add_subplot(len(band_list),1,1+m)
    big_M_0_trace = []
    for n in trace_run_nums:
        M_0_trace = loadtxt(data_dir + band_list[m] + "_" + plr_db_prefix + "_1" + str(n) + "_M_0_trace.txt")[-n_samples_analysis:]
                
        big_M_0_trace.append(M_0_trace.tolist())
        ax.plot(M_0_trace, color=colors_list[n-1], alpha=0.1)
    big_M_0_trace = array(big_M_0_trace)
    big_M_0_mean = big_M_0_trace.mean()
    big_M_0_std = big_M_0_trace.std()
    ax.hlines(big_M_0_mean, 0, n_samples_analysis, colors="k")
    ax.hlines([big_M_0_mean+big_M_0_std, big_M_0_mean-big_M_0_std], 0, n_samples_analysis, colors="k", linestyles="--")
    
    ylabel_obj = ax.set_ylabel(band_list[m] + "\n" + r" %.4f$\pm$%.4f" % (big_M_0_mean, big_M_0_std))
    ylabel_obj.set_rotation(0)
    ylabel_obj.set_horizontalalignment("right")
    ax.set_xlim(0, n_samples_analysis)
    # pos = [left, bottom, width, height]
    ax_pos = ax.get_position().bounds
    ax.set_position([ax_pos[0]+0.04, ax_pos[1]-0.0105*(m-3.8), ax_pos[2]+0.03, ax_pos[3]+0.015])
    
    if m < len(band_list)-1:
        for tl in ax.get_xticklabels():
            tl.set_visible(False)
fig.suptitle("M_0 traces for each band")
canvas = FigureCanvas(fig)
canvas.print_figure(plr_db_prefix + "_M_0_traces.png", dpi=300)
close("all")


fig = plt.figure(figsize=(8, 12))
for m in range(len(band_list)):
    ax = fig.add_subplot(len(band_list),1,1+m)
    big_alpha_trace = []
    for n in trace_run_nums:
        alpha_trace = loadtxt(data_dir + band_list[m] + "_" + plr_db_prefix + "_1" + str(n) + "_alpha_trace.txt")[-n_samples_analysis:]
                
        big_alpha_trace.append(alpha_trace.tolist())
        ax.plot(alpha_trace, color=colors_list[n-1], alpha=0.1)
    big_alpha_trace = array(big_alpha_trace)
    big_alpha_mean = big_alpha_trace.mean()
    big_alpha_std = big_alpha_trace.std()
    ax.hlines(big_alpha_mean, 0, n_samples_analysis, colors="k")
    ax.hlines([big_alpha_mean+big_alpha_std, big_alpha_mean-big_alpha_std], 0, n_samples_analysis, colors="k", linestyles="--")
    
    ylabel_obj = ax.set_ylabel(band_list[m] + "\n" + r" %.4f$\pm$%.4f" % (big_alpha_mean, big_alpha_std))
    ylabel_obj.set_rotation(0)
    ylabel_obj.set_horizontalalignment("right")
    ax.set_xlim(0, n_samples_analysis)
    # pos = [left, bottom, width, height]
    ax_pos = ax.get_position().bounds
    ax.set_position([ax_pos[0]+0.04, ax_pos[1]-0.0105*(m-3.8), ax_pos[2]+0.03, ax_pos[3]+0.015])
    
    if m < len(band_list)-1:
        for tl in ax.get_xticklabels():
            tl.set_visible(False)
fig.suptitle("alpha traces for each band")
canvas = FigureCanvas(fig)
canvas.print_figure(plr_db_prefix + "_alpha_traces.png", dpi=300)
close("all")