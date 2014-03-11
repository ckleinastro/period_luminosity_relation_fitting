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

band_wavelengths = {    "U":0.3663,
                        "B":0.4361,
                        "hipp":0.5170,
                        "V":0.5448,
                        "R":0.6407,
                        "I":0.7980,
                        "z":0.8896,
                        "J":1.22,
                        "H":1.63,
                        "K":2.19,
                        "W1":3.4,
                        "W2":4.6,
                        "W3":12.0
                    }


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

"""
0   name
1   type
2   blazhko
3   period
4   log10(P_f/P_0)
5   metallicity
6   prior_mu
7   prior_mu_err
8   SF_ebv
9   SF_ebv_err
10  m_U_obs
11  m_U_obs_err
12  m_B_obs
13  m_B_obs_err
14  m_hipp_obs
15  m_hipp_obs_err
16  m_V_obs
17  m_V_obs_err
18  m_R_obs
19  m_R_obs_err
20  m_I_obs
21  m_I_obs_err
22  m_z_obs
23  m_z_obs_err
24  m_J_obs
25  m_J_obs_err
26  m_H_obs
27  m_H_obs_err
28  m_K_obs
29  m_K_obs_err
30  m_W1_obs
31  m_W1_obs_err
32  m_W2_obs
33  m_W2_obs_err
34  m_W3_obs
35  m_W3_obs_err
36  post_mu
37  post_mu_err
38  post_ebv
39  post_ebv_err
40  M_U_fit
41  M_U_fit_err
42  M_B_fit
43  M_B_fit_err
44  M_hipp_fit
45  M_hipp_fit_err
46  M_V_fit
47  M_V_fit_err
48  M_R_fit
49  M_R_fit_err
50  M_I_fit
51  M_I_fit_err
52  M_z_fit
53  M_z_fit_err
54  M_J_fit
55  M_J_fit_err
56  M_H_fit
57  M_H_fit_err
58  M_K_fit
59  M_K_fit_err
60  M_W1_fit
61  M_W1_fit_err
62  M_W2_fit
63  M_W2_fit_err
64  M_W3_fit
65  M_W3_fit_err
"""


plr_fit_parameters_dtype=dtype([('band', str, 4),
                                ('M_0', 'float'),
                                ('M_0_err', 'float'),
                                ('alpha', 'float'),
                                ('alpha_err', 'float'),
                                ('sigma_intrinsic', 'float'),
                                ('sigma_intrinsic_err', 'float'),
                                ('sigma_instrumental', 'float'),
                                ('sigma_instrumental_err', 'float')])

plr_fit_parameters = loadtxt("plr_fit_parameters_table.txt", skiprows=1, dtype=plr_fit_parameters_dtype)



wavelengths = []
for band in band_list:
    wavelengths.append(band_wavelengths[band])
wavelengths = array(wavelengths)

log_wavelength_grid = linspace(log10(wavelengths)[0], log10(wavelengths)[-1], 1000)







spline_order = 2
spline_smoothness = 0.2

fig = plt.figure(figsize=(3.3, 3.3))
ax1 = subplot(111)




    
tck = interpolate.splrep(log10(wavelengths), plr_fit_parameters["alpha"], k=spline_order, s=spline_smoothness)
ax1.plot(log_wavelength_grid, interpolate.splev(log_wavelength_grid,tck,der=0), color="k", alpha=0.5)

spline_fit_slopes = interpolate.splev(log10(wavelengths),tck,der=0)

tck_up = interpolate.splrep(log10(wavelengths), spline_fit_slopes+plr_fit_parameters["alpha_err"], k=spline_order, s=0.06)
tck_down = interpolate.splrep(log10(wavelengths), spline_fit_slopes-plr_fit_parameters["alpha_err"], k=spline_order, s=0.06)
    
ax1.fill_between(log_wavelength_grid, interpolate.splev(log_wavelength_grid,tck_up,der=0), interpolate.splev(log_wavelength_grid,tck_down,der=0), alpha=0.25, color="k")

# Sollima 2006
ax1.errorbar(log10(band_wavelengths["K"]), -2.38, 0.04, color="limegreen", marker="s", markersize=3, linestyle="none", label="Sollima et al. (2006)")

# Catelan 2004
ax1.plot(log10(array([band_wavelengths["I"], band_wavelengths["J"], band_wavelengths["H"], band_wavelengths["K"]])), array([-1.132, -1.773, -2.313, -2.353]), color="blue", marker="h", markersize=3, linestyle="none", label="Catelan et al. (2004)")

# Madore 2013
ax1.errorbar(log10(array([band_wavelengths["W1"], band_wavelengths["W2"], band_wavelengths["W3"]])), array([-2.44, -2.55, -2.58]), array([0.95, 0.89, 0.97]), color="red", marker="D", markersize=3, linestyle="none", label="Madore et al. (2013)")

# Dambis 2014
ax1.errorbar(log10(array([band_wavelengths["W1"], band_wavelengths["W2"]])), array([-2.381, -2.269]), array([0.097, 0.127]), color="cyan", marker="v", markersize=3, linestyle="none", label="Dambis et al. (2014)")

ax1.errorbar(log10(wavelengths), plr_fit_parameters["alpha"], plr_fit_parameters["alpha_err"], color="k", marker="o", markersize=3, linestyle="none", label="This Paper")

ax1.legend(loc="lower right", numpoints=1, fontsize=10, handlelength=0.7, labelspacing=0.5, handletextpad=0.4)
ax1.set_xlim(log10(wavelengths).min()-0.1, log10(wavelengths).max()+0.1)
ax1.set_ylim(0.5, -2.9)

ax1.set_xlabel(r"$\log_{10}(\lambda [\mu{\rm m}])$")
ax1.set_ylabel(r"$\alpha$")



# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.2)
minorLocator_x = MultipleLocator(0.1)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.5)
minorLocator_y1 = MultipleLocator(0.25)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)


# pos =         [left, bottom, width, height]
ax1.set_position([0.17, 0.15, 0.81, 0.83])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "slope.pdf", dpi=300)
close("all")
