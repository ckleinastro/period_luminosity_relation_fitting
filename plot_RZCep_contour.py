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

P_0 = 0.52853966619770265

data_dir = "plots/"


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




full_EBV_trace = loadtxt("RZCep_EBV_trace.txt")
full_mu_trace = loadtxt("RZCep_mu_trace.txt")
full_RV_trace = loadtxt("RZCep_RV_trace.txt")


x_min = full_EBV_trace.mean()-5*full_EBV_trace.std()
x_max = full_EBV_trace.mean()+5*full_EBV_trace.std()
y_min_mu = full_mu_trace.mean()-5*full_mu_trace.std()
y_max_mu = full_mu_trace.mean()+5*full_mu_trace.std()

x_data = array([full_EBV_trace]).T
y_data_mu = array([full_mu_trace]).T
rvs_mu = np.append(x_data, y_data_mu, axis=1)

kde_mu = stats.kde.gaussian_kde(rvs_mu.T)

# Regular grid to evaluate kde_mu upon
x_flat = np.r_[x_min:x_max:256j]
y_flat_mu = np.r_[y_min_mu:y_max_mu:256j]
x,y_mu = np.meshgrid(x_flat,y_flat_mu)
grid_coords = np.append(x.reshape(-1,1),y_mu.reshape(-1,1),axis=1)

z_mu = kde_mu(grid_coords.T)
z_mu = z_mu.reshape(256,256)

l_999_mu = find_level(z_mu, 0.001)
l_99_mu = find_level(z_mu, 0.01)
l_975_mu = find_level(z_mu, 0.025)
l_95_mu = find_level(z_mu, 0.05)
l_90_mu = find_level(z_mu, 0.10)
l_85_mu = find_level(z_mu, 0.15)
l_80_mu = find_level(z_mu, 0.20)
l_70_mu = find_level(z_mu, 0.30)
levels_list_mu=[l_70_mu, l_80_mu, l_85_mu, l_90_mu, l_95_mu, l_975_mu, l_99_mu, l_999_mu, z_mu.max()]





y_min_RV = full_RV_trace.mean()-5*full_RV_trace.std()
y_max_RV = full_RV_trace.mean()+5*full_RV_trace.std()

y_data_RV = array([full_RV_trace]).T
rvs_RV = np.append(x_data, y_data_RV, axis=1)

kde_RV = stats.kde.gaussian_kde(rvs_RV.T)

# Regular grid to evaluate kde_RV upon
y_flat_RV = np.r_[y_min_RV:y_max_RV:256j]
x,y_RV = np.meshgrid(x_flat,y_flat_RV)
grid_coords = np.append(x.reshape(-1,1),y_RV.reshape(-1,1),axis=1)

z_RV = kde_RV(grid_coords.T)
z_RV = z_RV.reshape(256,256)



l_999_RV = find_level(z_RV, 0.001)
l_99_RV = find_level(z_RV, 0.01)
l_975_RV = find_level(z_RV, 0.025)
l_95_RV = find_level(z_RV, 0.05)
l_90_RV = find_level(z_RV, 0.10)
l_85_RV = find_level(z_RV, 0.15)
l_80_RV = find_level(z_RV, 0.20)
l_70_RV = find_level(z_RV, 0.30)
levels_list_RV=[l_70_RV, l_80_RV, l_85_RV, l_90_RV, l_95_RV, l_975_RV, l_99_RV, l_999_RV, z_RV.max()]







nullfmt = NullFormatter()   # no labels

# definitions for the axes
# left, width = 0.17, 0.58
# bottom, height = 0.17, 0.58
# bottom_h = left_h = left+width+0.02
# 
# rect_scatter = [left, bottom, width, height]
# rect_histx = [left, bottom_h, width, 0.2]
# rect_histy = [left_h, bottom, 0.2, height]
# 
# rect_band_label = [left_h, bottom_h, 0.2, 0.2]


fig_width = 3.3
fig_height = 5.15

scatter_size = 1.914

# pos =          [left, bottom, width, height]
# definitions for the axes
left, width = 0.17, scatter_size/fig_width
bottom, height_RV, height_mu = 0.10, scatter_size/fig_height, scatter_size/fig_height
left_h = left+width+0.066/fig_width
bottom_h = bottom+height_RV+height_mu + 0.066/fig_height

rect_scatter_mu = [left, bottom+height_RV, width, height_mu]
rect_scatter_RV = [left, bottom, width, height_RV]

rect_histx = [left, bottom_h, width, 0.66/fig_height]

rect_histy_mu = [left_h, bottom+height_RV, 0.66/fig_width, height_mu]
rect_histy_RV = [left_h, bottom, 0.66/fig_width, height_RV]

rect_band_label = [left_h, bottom_h, 0.66/fig_width, 0.66/fig_height]




fig = plt.figure(figsize=(fig_width, fig_height))

axScatter_mu = plt.axes(rect_scatter_mu)
axScatter_RV = plt.axes(rect_scatter_RV)
axHistx = plt.axes(rect_histx)
axHisty_mu = plt.axes(rect_histy_mu)
axHisty_RV = plt.axes(rect_histy_RV)

axBandLabel = plt.axes(rect_band_label)

axBandLabel.yaxis.set_ticks([])
axBandLabel.xaxis.set_ticks([])
axBandLabel.set_xlim(0,1)
axBandLabel.set_ylim(0,1)
axBandLabel.text(0.5, 0.5, "RZCep", fontsize=14, ha='center', va='center')



axScatter_mu.contourf(x, y_mu, z_mu, levels=levels_list_mu, origin="lower", cmap=cm.gist_yarg)
axScatter_mu.contour(x, y_mu, z_mu, levels=levels_list_mu, origin="lower", colors='k')
axScatter_mu.errorbar(x_data.mean(), y_data_mu.mean(), y_data_mu.std(), x_data.std(), marker="o", color="r")

axScatter_RV.contourf(x, y_RV, z_RV, levels=levels_list_RV, origin="lower", cmap=cm.gist_yarg)
axScatter_RV.contour(x, y_RV, z_RV, levels=levels_list_RV, origin="lower", colors='k')
axScatter_RV.errorbar(x_data.mean(), y_data_RV.mean(), y_data_RV.std(), x_data.std(), marker="o", color="r")

axScatter_RV.set_xlabel(r"$E(B-V)$")
axScatter_mu.set_ylabel(r"$\mu$")
axScatter_RV.set_ylabel(r"$R_V$")
axScatter_mu.set_xlim(x_min, x_max)
axScatter_RV.set_xlim(x_min, x_max)

axScatter_RV.set_ylim(y_min_RV, y_max_RV)
axScatter_mu.set_ylim(y_min_mu, y_max_mu)

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_y_mu = MultipleLocator(0.03)
minorLocator_y_mu = MultipleLocator(0.015)
axScatter_mu.yaxis.set_major_locator(majorLocator_y_mu)
axScatter_mu.yaxis.set_minor_locator(minorLocator_y_mu)

majorLocator_y_RV = MultipleLocator(1)
minorLocator_y_RV = MultipleLocator(0.5)
axScatter_RV.yaxis.set_major_locator(majorLocator_y_RV)
axScatter_RV.yaxis.set_minor_locator(minorLocator_y_RV)

majorLocator_x_EBV = MultipleLocator(0.1)
minorLocator_x_EBV = MultipleLocator(0.05)
axScatter_RV.xaxis.set_major_locator(majorLocator_x_EBV)
axScatter_RV.xaxis.set_minor_locator(minorLocator_x_EBV)




axHistx.hist(full_EBV_trace, bins=100, range=axScatter_RV.get_xlim(), normed=False, histtype="stepfilled", color="gray", alpha=1.0)

axHisty_mu.hist(full_mu_trace, bins=100, range=axScatter_mu.get_ylim(), normed=False, histtype="stepfilled", orientation="horizontal", color="gray", alpha=1.0)

axHisty_RV.hist(full_RV_trace, bins=100, range=axScatter_RV.get_ylim(), normed=False, histtype="stepfilled", orientation="horizontal", color="gray", alpha=1.0)

axHistx.set_ylabel("Count")
axHisty_RV.set_xlabel("Count")

axHistx.set_xlim( axScatter_RV.get_xlim() )
axHisty_mu.set_ylim( axScatter_mu.get_ylim() )
axHisty_RV.set_ylim( axScatter_RV.get_ylim() )

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty_RV.yaxis.set_major_formatter(nullfmt)
axHisty_mu.yaxis.set_major_formatter(nullfmt)
axScatter_mu.xaxis.set_major_formatter(nullfmt)
xticksHisty_mu = axHisty_mu.xaxis.get_major_ticks()
for x_tick_label in xticksHisty_mu:
    x_tick_label.label1.set_visible(False)

xticksHisty_RV = axHisty_RV.xaxis.get_major_ticks()
for x_tick_label in xticksHisty_RV:
    x_tick_label.label1.set_visible(False)

yticksHistx = axHistx.yaxis.get_major_ticks()
for y_tick_label in yticksHistx:
    y_tick_label.label1.set_visible(False)


canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "RZCep_prediction_contour_with_RV.pdf", dpi=300)
close("all")




