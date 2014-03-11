from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from numpy import array, log10, loadtxt
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





plr_data_dict = {}

for n in range(len(band_list)):
    band = band_list[n]
    linear_per = []
    fund_log_per = []
    M_band = []
    M_band_err = []
    type = []
    blazhko = []

    star_fit_data_file = file("full_rrl_fit_table.txt", "r")
    first_line = star_fit_data_file.readline()
    for line in star_fit_data_file:
        M = line.split()[40+2*n]
        if M != "---":
            linear_per.append(float(line.split()[3]))
            fund_log_per.append(float(line.split()[4]))
            M_band.append(float(M))
            M_band_err.append(float(line.split()[41+2*n]))
            type.append(line.split()[1])
            blazhko.append(line.split()[2]=="True")
        
        
    
    
    star_fit_data_file.close()
    
    linear_per = array(linear_per)
    fund_log_per = array(fund_log_per)
    M_band = array(M_band)
    M_band_err = array(M_band_err)
    type = array(type)
    blazhko = array(blazhko)

    plr_data_dict[band] = {}
    plr_data_dict[band]["type_array"] = type
    plr_data_dict[band]["blazhko_array"] = blazhko
    
    plr_data_dict[band]["linear_per_array"] = linear_per
    plr_data_dict[band]["fund_log_per_array"] = fund_log_per
    plr_data_dict[band]["M_band_array"] = M_band
    plr_data_dict[band]["M_band_err_array"] = M_band_err

    band_fit_parameters = plr_fit_parameters[plr_fit_parameters["band"] == band]
    
    plr_data_dict[band]["band_fit_parameters"] = band_fit_parameters



logper_grid_bounds = array([-0.55, -0.02])
normed_logper_grid_bounds = logper_grid_bounds - log10(P_0)
logper_grid = linspace(normed_logper_grid_bounds[0], normed_logper_grid_bounds[-1], 1000)

min_prediction_sigma_array = []

for band in band_list:

    M_0_fit = plr_data_dict[band]["band_fit_parameters"]["M_0"]
    M_0_err_fit = plr_data_dict[band]["band_fit_parameters"]["M_0_err"]
    alpha_fit = plr_data_dict[band]["band_fit_parameters"]["alpha"]
    alpha_err_fit = plr_data_dict[band]["band_fit_parameters"]["alpha_err"]

    full_M_0_trace = []
    full_alpha_trace = []
    for n in trace_run_nums:
        M_0_trace = loadtxt(data_dir + band + "_" + plr_db_prefix + "_1" + str(n) + "_M_0_trace.txt")
        alpha_trace = loadtxt(data_dir + band + "_" + plr_db_prefix + "_1" + str(n) + "_alpha_trace.txt")
        full_M_0_trace.append(M_0_trace.tolist())
        full_alpha_trace.append(alpha_trace.tolist())
        
    full_M_0_trace = array(full_M_0_trace)
    full_M_0_trace = full_M_0_trace.flatten()
    full_alpha_trace = array(full_alpha_trace)
    full_alpha_trace = full_alpha_trace.flatten()
    
    squared_m_errs = (plr_data_dict[band]["M_band_array"] - (M_0_fit + alpha_fit*plr_data_dict[band]["fund_log_per_array"]))**2
    sig_m = sqrt(squared_m_errs.mean())

    ci_err_grid = []
    for logper_val in logper_grid:
        ci_err_grid.append(sqrt(sig_m**2 + std(full_M_0_trace + full_alpha_trace * logper_val)**2))
    ci_err_grid = array(ci_err_grid)
    
    min_prediction_sigma = ci_err_grid.min()
    min_prediction_sigma_array.append(min_prediction_sigma)

min_prediction_sigma_array = array(min_prediction_sigma_array)


fig = plt.figure(figsize=(3.3, 4.5))
ax1 = subplot(111)

wavelengths = []
sigma_intrinsic_array = []
sigma_intrinsic_err_array = []
sigma_instrumental_array = []
sigma_instrumental_err_array = []
sigma_prediction_array = []

for m in range(len(band_list)):
    band = band_list[m]
    band_wavelength = band_wavelengths[band]
    sigma_intrinsic = plr_data_dict[band]["band_fit_parameters"]["sigma_intrinsic"][0]
    sigma_intrinsic_err = plr_data_dict[band]["band_fit_parameters"]["sigma_intrinsic_err"][0]
    sigma_instrumental = plr_data_dict[band]["band_fit_parameters"]["sigma_instrumental"][0]
    sigma_instrumental_err = plr_data_dict[band]["band_fit_parameters"]["sigma_instrumental_err"][0]
    
    sigma_prediction = min_prediction_sigma_array[m]
    
    wavelengths.append(band_wavelength)
    
    sigma_intrinsic_array.append(sigma_intrinsic)
    sigma_intrinsic_err_array.append(sigma_intrinsic_err)
    sigma_instrumental_array.append(sigma_instrumental)
    sigma_instrumental_err_array.append(sigma_instrumental_err)
    sigma_prediction_array.append(sigma_prediction)

    


sigma_intrinsic_array = array(sigma_intrinsic_array)
sigma_intrinsic_err_array = array(sigma_intrinsic_err_array)
sigma_instrumental_array = array(sigma_instrumental_array)
sigma_instrumental_err_array = array(sigma_instrumental_err_array)
sigma_prediction_array = array(sigma_prediction_array)
wavelengths = array(wavelengths)


# ax1.errorbar(log10(wavelengths), sigma_intrinsic_array, sigma_intrinsic_err_array, marker="o", linestyle="-", markersize=3, color="magenta", label="Intrinsic Scatter")
ax1.plot(log10(wavelengths), sigma_intrinsic_array, marker="o", linestyle="-", markersize=3, color="magenta", label="Intrinsic Scatter")
ax1.fill_between(log10(wavelengths), sigma_intrinsic_array+2*sigma_intrinsic_err_array, sigma_intrinsic_array-2*sigma_intrinsic_err_array, alpha=0.25, color="magenta")


# ax1.errorbar(log10(wavelengths), sigma_instrumental_array, sigma_instrumental_err_array, marker="o", linestyle="-", markersize=3, color="green", label="Photometric Error")
ax1.plot(log10(wavelengths), sigma_instrumental_array, marker="o", linestyle="-", markersize=3, color="green", label="Photometric Error")
ax1.fill_between(log10(wavelengths), sigma_instrumental_array+2*sigma_instrumental_err_array, sigma_instrumental_array-2*sigma_instrumental_err_array, alpha=0.25, color="green")

# ax1.plot(log10(wavelengths), sigma_prediction_array, marker="o", linestyle="-", markersize=3, color="orange", label="Prediction Uncertainty")

ax1.text(0.07, 0.087, r"Shaded regions are $\pm 2\sigma$.", fontsize=10)

ax1.legend(loc="upper right", numpoints=1, fontsize=10, handlelength=0.7, labelspacing=0.5, handletextpad=0.4)
ax1.set_xlim(log10(band_wavelengths["U"])-0.1, log10(band_wavelengths["W3"])+0.1)

ax1.set_ylim(0, 0.119)


# ax1.set_xlabel(r"$\lambda^{-1}$ [$\mu{\rm m}^{-1}$]")
ax1.set_xlabel(r"$\log_{10}(\lambda [\mu{\rm m}])$")
ax1.set_ylabel(r"Magnitudes")



# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.2)
minorLocator_x = MultipleLocator(0.1)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.02)
minorLocator_y1 = MultipleLocator(0.01)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)


# pos =         [left, bottom, width, height]
ax1.set_position([0.18, 0.12, 0.80, 0.87])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "sigmas.pdf", dpi=300)
close("all")
