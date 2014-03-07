import pymc
from scipy import loadtxt, log10, ones, sqrt, mean, polyfit, linspace, random
from scipy import zeros, append, insert, hstack, vstack, concatenate, savetxt
from scipy import identity, dot
from scipy import interpolate
from numpy import linalg, sin, cos, arcsin, savetxt
import sys
from pylab import *
import operator
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import cPickle as pickle

import matplotlib.cm as cm
import matplotlib.colors as mcolors

if len(sys.argv) > 1:
    plr_db_num = sys.argv[1]
else:
    plr_db_num = "test"
    

const_R_V = 3.1

# to reload a database:
import socket
if socket.gethostname() == "Christopher-Kleins-MacBook-Pro.local":
    M = pymc.database.pickle.load('/Volumes/Extra_HDD/Research/Multiband_PLR_Fitting/plr_' + plr_db_num + '.pickle')
    plot_dir = "/Volumes/Extra_HDD/Research/Multiband_PLR_Fitting/trace_plots/"
    # data_dir = "/Volumes/Extra_HDD/Research/Multiband_PLR_Fitting/trace_data/"
    data_dir = "/Volumes/Extra_HDD/Research/Multiband_PLR_Fitting/lowexpandedEBV_trace_data/"
    # plot_dir = "trace_plots/"
elif socket.gethostname() == "anathem":
    M = pymc.database.pickle.load('/big_data/cklein/plr_fitting_traces/plr_' + plr_db_num + '.pickle')
    plot_dir = "/big_data/cklein/plr_fitting_traces/trace_plots/"
    data_dir = "/big_data/cklein/plr_fitting_traces/trace_data/"
else:
    print "Error, not running on known host, exiting."
    sys.exit()


n_samples = M.M_0.gettrace().shape[0]
n_samples_analysis = int(n_samples/2.0)


rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=18)

band_list = ["U", "B", "hipp", "V", "R", "I", "z", "J", "H", "K", "W1", "W2", "W3"]
# band_list = ["W1", "W2"]

test_names = array(['AACMi', 'ABUMa', 'AEBoo', 'AFVel', 'AFVir', 'AMTuc', 'AMVir',
       'ANSer', 'APSer', 'ARHer', 'ATAnd', 'ATVir',
       'AUVir', 'AVPeg', 'AVVir', 'AXLeo', 'BBEri', 'BCDra', 'BHPeg',
       'BKDra', 'BNPav', 'BPPav', 'BRAqr', 'BTDra', 'BVAqr',
       'BXLeo', 'CGLib', 'CGPeg', 'CIAnd', 'CSEri', 'DDHya', 'DHPeg',
       'DNAqr', 'DXDel', 'FWLup', 'HHPup', 'HKPup', 'IKHya', 'IOLyr',
       'MSAra', 'MTTel', 'RRCet', 'RRGem', 'RRLeo', 'RRLyr', 'RSBoo',
       'RUCet', 'RUPsc', 'RUScl', 'RVCap', 'RVCet', 'RVCrB', 'RVOct',
       'RVUMa', 'RWCnc', 'RWDra', 'RXCet', 'RXCol', 'RXEri', 'RYCol',
       'RYOct', 'RZCet', 'RZCVn', 'SAra', 'SCom', 'SSCVn',
       'SSFor', 'SSLeo', 'SSOct', 'STBoo', 'STCom', 'STCVn', 'STLeo',
       'STVir', 'SUDra', 'SVEri', 'SVHya', 'SVScl', 'SWAnd', 'SWAqr',
       'SWDra', 'SXAqr', 'SXFor', 'SXUMa', 'SZGem', 'TSex', 'TTCnc',
       'TTLyn', 'TUUMa', 'TVBoo', 'TVCrB', 'TWBoo', 'TWHer', 'TWLyn',
       'TYAps', 'TZAur', 'UCom', 'ULep', 'UPic', 'UUCet', 'UUVir', 'UVOct',
       'UYBoo', 'UYCyg', 'UZCVn', 'V341Aql', 'V413CrA', 'V440Sgr',
       'V445Oph', 'V499Cen', 'V675Sgr', 'VInd', 'VWScl', 'VXHer', 'VXScl',
       'VYLib', 'VYSer', 'VZHer', 'VZPeg', 'WCrt', 'WCVn', 'WTuc', 'WYAnt',
       'WYPav', 'WZHya', 'XAri', 'XCrt', 'XXAnd', 'XXPup', 'XZAps',
       'XZCyg', 'XZDra', 'YZCap', 'ZMic'], 
      dtype='|S7')

all_names = test_names





rrl_info_type=dtype([('name', str, 7),
                 ('type', str, 4),
                 ('period', 'float'),
                 ('metallicity', 'float')])
rrl_info = loadtxt("rrl_info.txt", dtype=rrl_info_type)
# all_names = rrl_info["name"]
short_rrl_info = []
for n in range(len(rrl_info)):
    if rrl_info["name"][n] in all_names:
        short_rrl_info.append(rrl_info[n])
rrl_info = array(short_rrl_info, dtype=rrl_info_type)



prior_mu_dtype=dtype([('name', str, 7),
                 ('mu', 'float'),
                 ('mu_err', 'float')])
prior_mu_info = loadtxt("prior_distance_moduli.txt", dtype=prior_mu_dtype)
short_prior_mu_info = []
for n in range(len(prior_mu_info)):
    if prior_mu_info["name"][n] in all_names:
        short_prior_mu_info.append(prior_mu_info[n])
prior_mu_info = array(short_prior_mu_info, dtype=prior_mu_dtype)




# Replace the metallacity-based mu prior with HST parallax where available
parallax_data = loadtxt("benedict_parallax.txt", skiprows=1, dtype=str)
for n in range(len(prior_mu_info['name'])):
    for hst_parallax_entry in parallax_data:
        if prior_mu_info['name'][n] == hst_parallax_entry[0]:
            prior_mu_info['mu'][n] = hst_parallax_entry[3]
            prior_mu_info['mu_err'][n] = hst_parallax_entry[4]


ebv_color_excess = []
# name, ra, dec, E(B-V), err_E(B-V)
sfd_extinction_file = file("extinction.tbl.txt")
for line in sfd_extinction_file:
    if line[0] not in ("|", "\\"):
        ebv_color_excess.append((line.split()[0], 
                                 float(line.split()[1]), float(line.split()[2]), 
                                 float(line.split()[5]), float(line.split()[6])))
sfd_extinction_file.close()
extinction_dtype = dtype([('name', str, 7),
                          ('ra', 'float'),
                          ('dec', 'float'),
                          ('ebv', 'float'),
                          ('ebv_err', 'float')])
ebv_color_excess = array(ebv_color_excess, dtype=extinction_dtype)

blazhko_stars = loadtxt("blazhko_list.txt", dtype=str)

mags = pickle.load( open( "mean_flux_mags_dict.p", "rb" ) )

observed_apparent_mags = []
observed_apparent_mag_errs = []
for band in band_list:
    for name in all_names:
        if band in mags[name].keys():
            observed_apparent_mags.append(mags[name][band][0])
            observed_apparent_mag_errs.append(mags[name][band][1])
observed_apparent_mags = array(observed_apparent_mags)
observed_apparent_mag_errs = array(observed_apparent_mag_errs)

fundamentalized_log_periods = append(log10(rrl_info["period"][rrl_info["type"]=="RRc"]) + 0.127, 
    log10(rrl_info["period"][rrl_info["type"]=="RRab"]))
P_0 = (10**fundamentalized_log_periods).mean()
for n in range(len(rrl_info["period"])):
    if rrl_info["type"][n] == "RRab":
        rrl_info["period"][n] = log10(rrl_info["period"][n]) - log10(P_0)
    elif rrl_info["type"][n] == "RRc":
        rrl_info["period"][n] = log10(rrl_info["period"][n]) + 0.127 - log10(P_0)
    else:
        print "Screwed up period transform for", rrl_info[n]



ebv_values = []
ebv_errs = []
for name in all_names:
    ebv_values.append(ebv_color_excess["ebv"][ebv_color_excess["name"]==name][0])
    ebv_errs.append(ebv_color_excess["ebv_err"][ebv_color_excess["name"]==name][0])
ebv_values = array(ebv_values)
ebv_errs = array(ebv_errs)




def ccm_extinction_law(wavelength):
    # Extinction law (A_lambda/A_V) from Cardelli, Clayton, and Mathis 1989
    x = 1.0/wavelength
    if x <= 1.1:    # If infrared or longer wavelength
        a = 0.574*x**1.61
        b = -0.527*x**1.61
    if x > 1.1:     # If optical wavelength 
        y = x-1.82
        a = 1.0 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    return a, b

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

extinction_ccm_a = []
extinction_ccm_b = []
for band in band_list:
    extinction_coefficients = ccm_extinction_law(band_wavelengths[band])
    extinction_ccm_a.append(extinction_coefficients[0])
    extinction_ccm_b.append(extinction_coefficients[1])







# sigma_post = M.sigma.gettrace()[-n_samples:]

# sigma_by_band_post = M.sigma_by_band.gettrace()[-n_samples:][:,10]
# 
# W1_M_0 = M.M_0.gettrace()[-n_samples:][:,10]
# 
# W1_alpha = M.alpha.gettrace()[-n_samples:][:,10]
# 
# EBV_post = M.EBV.gettrace()[-n_samples:][:,0]

"""
The adaptive scale factors are tuning parameters for the Metropolis step method.
They should always be diminishing towards zero in their traces. After the burn-in
period, the tuning parameters should change infrequently and slightly (if at all).
"""
# plot(M.Metropolis_sigma_by_band_adaptive_scale_factor.gettrace())
# plot(M.Metropolis_M_0_adaptive_scale_factor.gettrace())
# plot(M.Metropolis_alpha_adaptive_scale_factor.gettrace())
# plot(M.Metropolis_EBV_adaptive_scale_factor.gettrace())
# plot(M.Metropolis_mus_adaptive_scale_factor.gettrace())
# plot(M.deviance.gettrace()[-n_samples_analysis:])

output_plr_file = file(plr_db_num + "_plr_fit_results.txt", "w")

output_plr_file.write("#band\tM_0\t\tM_0_err\talpha\talpha_err\tsigma\tsigma_err\n")

plr_fitted_values = []

for n in range(len(band_list)):
    if band_list[n] != "hipp":
        tab_str = "\t\t"
    else:
        tab_str = "\t"
    output_plr_file.write("%s%s%.4f\t%.4f\t%.4f\t%.4f\t\t%.4f\t%.4f\n" % 
        (   band_list[n], tab_str, 
            M.M_0.gettrace()[-n_samples_analysis:][:,n].mean(),
            M.M_0.gettrace()[-n_samples_analysis:][:,n].std(),
            M.alpha.gettrace()[-n_samples_analysis:][:,n].mean(),
            M.alpha.gettrace()[-n_samples_analysis:][:,n].std(),
            M.sigma_by_band.gettrace()[-n_samples_analysis:][:,n].mean(),
            M.sigma_by_band.gettrace()[-n_samples_analysis:][:,n].std() ))
    
    plr_fitted_values.append([ band_list[n], 
        band_wavelengths[band_list[n]],
        M.M_0.gettrace()[-n_samples_analysis:][:,n].mean(),
        M.M_0.gettrace()[-n_samples_analysis:][:,n].std(),
        M.alpha.gettrace()[-n_samples_analysis:][:,n].mean(),
        M.alpha.gettrace()[-n_samples_analysis:][:,n].std()
        ])
    
    savetxt(data_dir + band_list[n] + "_" + plr_db_num + "_M_0_trace.txt", M.M_0.gettrace()[-n_samples_analysis:][:,n], fmt='%.7f')
    savetxt(data_dir + band_list[n] + "_" + plr_db_num + "_alpha_trace.txt", M.alpha.gettrace()[-n_samples_analysis:][:,n], fmt='%.7f')
    savetxt(data_dir + band_list[n] + "_" + plr_db_num + "_sigma_trace.txt", M.sigma_by_band.gettrace()[-n_samples_analysis:][:,n], fmt='%.7f')
    
output_plr_file.close()





for n in range(len(all_names[:])):
    savetxt(data_dir + all_names[n] + "_" + plr_db_num + "_mu_trace.txt", M.mus.gettrace()[-n_samples_analysis:][:,n], fmt='%.7f')
    savetxt(data_dir + all_names[n] + "_" + plr_db_num + "_ebv_trace.txt", M.EBV.gettrace()[-n_samples_analysis:][:,n], fmt='%.7f')










sys.exit()






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

color_vals = linspace(1, 0, len(band_list))
colors_list = []
for cv in color_vals:
    colors_list.append(make_color(cv))

# for n in range(len(band_list)):
#     scatter(M.M_0.gettrace()[-n_samples_analysis:][:,n], M.alpha.gettrace()[-n_samples_analysis:][:,n], color=colors_list[n], alpha=0.1)
# 

"""
U		1.0602	0.0586	-0.4927	0.6497		0.2421	0.0437
B		0.8222	0.0251	-0.2831	0.2991		0.0989	0.0168
hipp	0.6551	0.0133	-0.5483	0.1633		0.0659	0.0071
V		0.5176	0.0182	-0.6227	0.2144		0.0646	0.0108
R		0.3370	0.0143	-0.9359	0.1748		0.0458	0.0090
I		0.1470	0.0291	-1.1259	0.3131		0.0610	0.0196
z		0.5965	0.0502	-1.3195	0.6553		0.1077	0.0385
J		-0.1306	0.0120	-1.8459	0.1575		0.0401	0.0079
H		-0.3396	0.0111	-2.3348	0.1535		0.0306	0.0071
K		-0.3409	0.0119	-2.6100	0.1557		0.0487	0.0084
W1		-0.4697	0.0065	-2.3359	0.0855		0.0036	0.0020
W2		-0.4597	0.0064	-2.3741	0.0854		0.0050	0.0019
W3		-0.4962	0.0073	-2.4372	0.1001		0.0226	0.0037
"""

approx_sigmas = [0.25, 0.1, 0.07, 0.07, 0.05, 0.07, 0.11, 0.05, 0.04, 0.05, 0.004, 0.0050, 0.03]
approx_M_0s = [1.0602, 0.8222, 0.6551, 0.5176, 0.3370, 0.1470, 0.5965, -0.1306, -0.3396, -0.3409, -0.4697, -0.4597, -0.4962]
approx_M_0s_widths = [0.0586, 0.0251, 0.0133, 0.0182, 0.0143, 0.0291, 0.0502, 0.0120, 0.0111, 0.0119, 0.0065, 0.0064, 0.0073]
approx_alphas = [-0.4927, -0.2831, -0.40, -0.6227, -0.9359, -1.1259, -1.3195, -1.8459, -2.3348, -2.6100, -2.3359, -2.3741, -2.4372]
approx_alphas_widths = [0.6497, 0.2991, 0.1633, 0.2144, 0.1748, 0.3131, 0.6553, 0.1575, 0.1535, 0.1557, 0.0855, 0.0854, 0.1001]

for m in range(len(band_list)):
    fig = plt.figure(figsize=(8, 12))
    ax1 = fig.add_subplot(6,1,1)
    ax1.plot(arange(n_samples-n_samples_analysis), M.sigma_by_band.gettrace()[:n_samples-n_samples_analysis][:,m], color="blue")
    ax1.plot(arange(n_samples-n_samples_analysis, n_samples), M.sigma_by_band.gettrace()[n_samples-n_samples_analysis:n_samples][:,m], color="red")
    ax2 = fig.add_subplot(6,1,2)
    ax2.plot(arange(n_samples-n_samples_analysis, n_samples), M.sigma_by_band.gettrace()[n_samples-n_samples_analysis:n_samples][:,m], color="red")
    sigma_mean_val = M.sigma_by_band.gettrace()[n_samples-n_samples_analysis:n_samples][:,m].mean()
    sigma_std_val = M.sigma_by_band.gettrace()[n_samples-n_samples_analysis:n_samples][:,m].std()
    
    ax2.plot([n_samples-n_samples_analysis, n_samples], [sigma_mean_val, sigma_mean_val], color="purple")
    ax2.plot([n_samples-n_samples_analysis, n_samples], [sigma_mean_val-sigma_std_val, sigma_mean_val-sigma_std_val], color="purple", linestyle="--")
    ax2.plot([n_samples-n_samples_analysis, n_samples], [sigma_mean_val+sigma_std_val, sigma_mean_val+sigma_std_val], color="purple", linestyle="--")
    
    ax1.set_xlim(0, n_samples)
    ax2.set_xlim(n_samples-n_samples_analysis, n_samples)
    ax1.set_ylim(0, 0.5)
    ax2.set_ylim(0, approx_sigmas[m]*2)
    
    ax1.set_title(band_list[m] + r" $\sigma$ Trace   analysis mean = %.4f +/- %.4f" % (sigma_mean_val, sigma_std_val))
    ax1.set_ylabel(r"$\sigma$")
    ax2.set_ylabel(r"$\sigma$ (zoomed in)")


    ax3 = fig.add_subplot(6,1,3)
    ax3.plot(arange(n_samples-n_samples_analysis), M.M_0.gettrace()[:n_samples-n_samples_analysis][:,m], color="blue")
    ax3.plot(arange(n_samples-n_samples_analysis, n_samples), M.M_0.gettrace()[n_samples-n_samples_analysis:n_samples][:,m], color="red")
    ax4 = fig.add_subplot(6,1,4)
    ax4.plot(arange(n_samples-n_samples_analysis, n_samples), M.M_0.gettrace()[n_samples-n_samples_analysis:n_samples][:,m], color="red")
    M_0_mean_val = M.M_0.gettrace()[n_samples-n_samples_analysis:n_samples][:,m].mean()
    M_0_std_val = M.M_0.gettrace()[n_samples-n_samples_analysis:n_samples][:,m].std()
    
    ax4.plot([n_samples-n_samples_analysis, n_samples], [M_0_mean_val, M_0_mean_val], color="purple")
    ax4.plot([n_samples-n_samples_analysis, n_samples], [M_0_mean_val-M_0_std_val, M_0_mean_val-M_0_std_val], color="purple", linestyle="--")
    ax4.plot([n_samples-n_samples_analysis, n_samples], [M_0_mean_val+M_0_std_val, M_0_mean_val+M_0_std_val], color="purple", linestyle="--")
    
    ax3.set_xlim(0, n_samples)
    ax4.set_xlim(n_samples-n_samples_analysis, n_samples)
    ax3.set_ylim(approx_M_0s[m]-6*approx_M_0s_widths[m], approx_M_0s[m]+6*approx_M_0s_widths[m])
    ax4.set_ylim(approx_M_0s[m]-3*approx_M_0s_widths[m], approx_M_0s[m]+3*approx_M_0s_widths[m])
    ax3.set_title(band_list[m] + r" $M_0$ Trace   analysis mean = %.4f +/- %.4f" % (M_0_mean_val, M_0_std_val))
    ax3.set_ylabel(r"$M_0$")
    ax4.set_ylabel(r"$M_0$ (zoomed in)")
    
    
    ax5 = fig.add_subplot(6,1,5)
    ax5.plot(arange(n_samples-n_samples_analysis), M.alpha.gettrace()[:n_samples-n_samples_analysis][:,m], color="blue")
    ax5.plot(arange(n_samples-n_samples_analysis, n_samples), M.alpha.gettrace()[n_samples-n_samples_analysis:n_samples][:,m], color="red")
    ax6 = fig.add_subplot(6,1,6)
    ax6.plot(arange(n_samples-n_samples_analysis, n_samples), M.alpha.gettrace()[n_samples-n_samples_analysis:n_samples][:,m], color="red")
    alpha_mean_val = M.alpha.gettrace()[n_samples-n_samples_analysis:n_samples][:,m].mean()
    alpha_std_val = M.alpha.gettrace()[n_samples-n_samples_analysis:n_samples][:,m].std()
    
    ax6.plot([n_samples-n_samples_analysis, n_samples], [alpha_mean_val, alpha_mean_val], color="purple")
    ax6.plot([n_samples-n_samples_analysis, n_samples], [alpha_mean_val-alpha_std_val, alpha_mean_val-alpha_std_val], color="purple", linestyle="--")
    ax6.plot([n_samples-n_samples_analysis, n_samples], [alpha_mean_val+alpha_std_val, alpha_mean_val+alpha_std_val], color="purple", linestyle="--")
    
    ax5.set_xlim(0, n_samples)
    ax6.set_xlim(n_samples-n_samples_analysis, n_samples)
    ax5.set_ylim(approx_alphas[m]-6*approx_alphas_widths[m], approx_alphas[m]+6*approx_alphas_widths[m])
    ax6.set_ylim(approx_alphas[m]-3*approx_alphas_widths[m], approx_alphas[m]+3*approx_alphas_widths[m])
    ax5.set_title(band_list[m] + r" $\alpha$ Trace   analysis mean = %.4f +/- %.4f" % (alpha_mean_val, alpha_std_val))
    ax6.set_xlabel("Trace Iteration")
    ax5.set_ylabel(r"$\alpha$")
    ax6.set_ylabel(r"$\alpha$ (zoomed in)")
    
    panel_height = 0.10666
    panel_width = 0.83
    left_pad = 0.12
    bottom_pad = 0.06
    vspacing_pad = 0.030
    
    # pos =          [left, bottom, width, height]
    ax1.set_position([left_pad, bottom_pad+5*(panel_height+vspacing_pad*1)+4*vspacing_pad,   panel_width,  panel_height])
    ax2.set_position([left_pad, bottom_pad+4*(panel_height+vspacing_pad*2),   panel_width,  panel_height])
    
    ax3.set_position([left_pad, bottom_pad+3*(panel_height+vspacing_pad*1)+2*vspacing_pad,   panel_width,  panel_height])
    ax4.set_position([left_pad, bottom_pad+2*(panel_height+vspacing_pad*2),   panel_width,  panel_height])
    
    ax5.set_position([left_pad, bottom_pad+1*(panel_height+vspacing_pad),   panel_width,  panel_height])
    ax6.set_position([left_pad, bottom_pad,   panel_width,  panel_height])
    
    canvas = FigureCanvas(fig)
    canvas.print_figure(plot_dir + band_list[m] + "_" + plr_db_num + "_traces.png", dpi=300)
    close("all")





# Only plot the first 12 star traces
for n in range(len(all_names[:])):
    fig = plt.figure(figsize=(8, 12*2/3.))
    ax1 = fig.add_subplot(4,1,1)
    ax1.plot(arange(n_samples-n_samples_analysis), M.mus.gettrace()[:n_samples-n_samples_analysis][:,n], color="blue")
    ax1.plot(arange(n_samples-n_samples_analysis, n_samples), M.mus.gettrace()[n_samples-n_samples_analysis:n_samples][:,n], color="red")
    ax2 = fig.add_subplot(4,1,2)
    ax2.plot(arange(n_samples-n_samples_analysis, n_samples), M.mus.gettrace()[n_samples-n_samples_analysis:n_samples][:,n], color="red")
    mu_mean_val = M.mus.gettrace()[n_samples-n_samples_analysis:n_samples][:,n].mean()
    mu_std_val = M.mus.gettrace()[n_samples-n_samples_analysis:n_samples][:,n].std()
    
    prior_mu = prior_mu_info[n][1]
    prior_mu_err = prior_mu_info[n][2]
    ax1.plot([0, n_samples], [prior_mu, prior_mu], color="green")
    ax1.plot([0, n_samples], [prior_mu+prior_mu_err, prior_mu+prior_mu_err], color="green",  linestyle="--")
    ax1.plot([0, n_samples], [prior_mu-prior_mu_err, prior_mu-prior_mu_err], color="green",  linestyle="--")
    
    ax2.plot([n_samples-n_samples_analysis, n_samples], [mu_mean_val, mu_mean_val], color="purple")
    ax2.plot([n_samples-n_samples_analysis, n_samples], [mu_mean_val-mu_std_val, mu_mean_val-mu_std_val], color="purple", linestyle="--")
    ax2.plot([n_samples-n_samples_analysis, n_samples], [mu_mean_val+mu_std_val, mu_mean_val+mu_std_val], color="purple", linestyle="--")
    
    ax1.set_xlim(0, n_samples)
    ax2.set_xlim(n_samples-n_samples_analysis, n_samples)
    
    ax1.set_ylim(prior_mu+2*prior_mu_err, prior_mu-2*prior_mu_err)
    ax2.set_ylim(mu_mean_val+3*mu_std_val, mu_mean_val-3*mu_std_val)
    
    ax1.set_title(all_names[n] + r" $\mu$ Trace   analysis mean = %.4f +/- %.4f" % (mu_mean_val, mu_std_val))
    ax1.set_ylabel(r"$\mu$")
    ax2.set_ylabel(r"$\mu$ (zoomed in)")


    ax3 = fig.add_subplot(4,1,3)
    ax3.plot(arange(n_samples-n_samples_analysis), M.EBV.gettrace()[:n_samples-n_samples_analysis][:,n], color="blue")
    ax3.plot(arange(n_samples-n_samples_analysis, n_samples), M.EBV.gettrace()[n_samples-n_samples_analysis:n_samples][:,n], color="red")
    ax4 = fig.add_subplot(4,1,4)
    ax4.plot(arange(n_samples-n_samples_analysis, n_samples), M.EBV.gettrace()[n_samples-n_samples_analysis:n_samples][:,n], color="red")
    EBV_mean_val = M.EBV.gettrace()[n_samples-n_samples_analysis:n_samples][:,n].mean()
    EBV_std_val = M.EBV.gettrace()[n_samples-n_samples_analysis:n_samples][:,n].std()
    
    ax4.plot([n_samples-n_samples_analysis, n_samples], [EBV_mean_val, EBV_mean_val], color="purple")
    ax4.plot([n_samples-n_samples_analysis, n_samples], [EBV_mean_val-EBV_std_val, EBV_mean_val-EBV_std_val], color="purple", linestyle="--")
    ax4.plot([n_samples-n_samples_analysis, n_samples], [EBV_mean_val+EBV_std_val, EBV_mean_val+EBV_std_val], color="purple", linestyle="--")
    
    sfd_ebv = ebv_values[n]
    
    
    ax3.set_xlim(0, n_samples)
    ax4.set_xlim(n_samples-n_samples_analysis, n_samples)
    
    ax3.set_ylim(0, max(2.6*sfd_ebv, 0.13, M.EBV.gettrace()[:][:,n].max()))
    
    ax3.set_title(all_names[n] + r" $E(B-V)$ Trace   analysis mean = %.4f +/- %.4f" % (EBV_mean_val, EBV_std_val))
    ax4.set_xlabel("Trace Iteration")
    ax3.set_ylabel(r"$E(B-V)$")
    ax4.set_ylabel(r"$E(B-V)$ (zoomed in)")
    
    panel_height = 0.185
    panel_width = 0.82
    left_pad = 0.13
    bottom_pad = 0.065
    vspacing_pad = 0.030
    
    # pos =          [left, bottom, width, height]
    
    ax1.set_position([left_pad, bottom_pad+3*(panel_height+vspacing_pad*1)+2*vspacing_pad,   panel_width,  panel_height])
    ax2.set_position([left_pad, bottom_pad+2*(panel_height+vspacing_pad*2),   panel_width,  panel_height])
    
    ax3.set_position([left_pad, bottom_pad+1*(panel_height+vspacing_pad),   panel_width,  panel_height])
    ax4.set_position([left_pad, bottom_pad,   panel_width,  panel_height])
    
    canvas = FigureCanvas(fig)
    canvas.print_figure(plot_dir + "star_traces/" +  all_names[n] + "_" + plr_db_num + "_traces.png", dpi=300)
    close("all")


def plot_tuning(parameter_name, trace):
    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(arange(2000, n_samples-n_samples_analysis), trace[2000:n_samples-n_samples_analysis], color="blue")
    ax1.plot(arange(n_samples-n_samples_analysis, n_samples), trace[n_samples-n_samples_analysis:n_samples], color="red")
    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(arange(n_samples-n_samples_analysis, n_samples), trace[n_samples-n_samples_analysis:n_samples], color="red")
    deviance_mean_val = trace[n_samples-n_samples_analysis:n_samples].mean()
    deviance_std_val = trace[n_samples-n_samples_analysis:n_samples].std()

    ax2.plot([n_samples-n_samples_analysis, n_samples], [deviance_mean_val, deviance_mean_val], color="purple")
    ax2.plot([n_samples-n_samples_analysis, n_samples], [deviance_mean_val-deviance_std_val, deviance_mean_val-deviance_std_val], color="purple", linestyle="--")
    ax2.plot([n_samples-n_samples_analysis, n_samples], [deviance_mean_val+deviance_std_val, deviance_mean_val+deviance_std_val], color="purple", linestyle="--")
    ax1.set_xlim(2000, n_samples)
    ax2.set_xlim(n_samples-n_samples_analysis, n_samples)
    ax1.set_title(parameter_name + " Trace   analysis mean = %.4f +/- %.4f" % (deviance_mean_val, deviance_std_val))
    ax1.set_ylabel(parameter_name)
    ax2.set_ylabel(parameter_name + " (zoomed in)")

    canvas = FigureCanvas(fig)
    canvas.print_figure(plot_dir + parameter_name + "_" + plr_db_num + "_trace.png", dpi=300)
    close("all")


plot_tuning("deviance", M.deviance.gettrace())
plot_tuning("sigma_SF", M.Metropolis_sigma_by_band_adaptive_scale_factor.gettrace())
plot_tuning("alpha_SF", M.Metropolis_alpha_adaptive_scale_factor.gettrace())
plot_tuning("EBV_SF", M.Metropolis_EBV_adaptive_scale_factor.gettrace())
plot_tuning("M_0_SF", M.Metropolis_M_0_adaptive_scale_factor.gettrace())
plot_tuning("mus_SF", M.Metropolis_mus_adaptive_scale_factor.gettrace())

    
post_ebv_values = []
post_ebv_errs = []
# post_RV_values = []
# post_RV_errs = []
post_mu_values = []
post_mu_errs = []

output_stars_file = file(plr_db_num + "_calibrated_rrls.txt", "w")
output_stars_file.write("#P_0 = %.4f\n" % P_0)
# output_stars_file.write("#name\tLog[P/P_0]\ttype\tblazhko\tR_V_post\tR_V_post_err\tebv_prior\tebv_prior_err\tebv_post\tebv_post_err\tmu_prior\tmu_prior_err\tmu_post\tmu_post_err\n")
output_stars_file.write("#name\tLog[P/P_0]\ttype\tblazhko\tebv_prior\tebv_prior_err\tebv_post\tebv_post_err\tmu_prior\tmu_prior_err\tmu_post\tmu_post_err\n")
for n in range(len(all_names)):
    blazhko_flag = "False"
    if rrl_info["name"][n] in blazhko_stars:
        blazhko_flag = "True"
    
    # post_RV_values.append(M.R_V.gettrace()[-n_samples_analysis:][:,n].mean())
    # post_RV_errs.append(M.R_V.gettrace()[-n_samples_analysis:][:,n].std())
    post_ebv_values.append(M.EBV.gettrace()[-n_samples_analysis:][:,n].mean())
    post_ebv_errs.append(M.EBV.gettrace()[-n_samples_analysis:][:,n].std())
    post_mu_values.append(M.mus.gettrace()[-n_samples_analysis:][:,n].mean())
    post_mu_errs.append(M.mus.gettrace()[-n_samples_analysis:][:,n].std())
    type_tab = "\t"
    if rrl_info["type"][n] == "RRc": type_tab = "\t\t"
    # output_stars_file.write("%s\t%.4f\t\t%s%s%s\t%.4f\t\t%.4f\t\t\t%.4f\t\t%.4f\t\t\t%.4f\t\t%.4f\t\t\t%.4f\t\t%.4f\t\t\t%.4f\t%.4f\n" %
    output_stars_file.write("%s\t%.4f\t\t%s%s%s\t%.4f\t\t%.4f\t\t\t%.4f\t\t%.4f\t\t\t%.4f\t\t%.4f\t\t\t%.4f\t%.4f\n" %  
        (   rrl_info["name"][n],
            rrl_info["period"][n],
            rrl_info["type"][n], type_tab,
            blazhko_flag,
            # M.R_V.gettrace()[-n_samples_analysis:][:,n].mean(),
            # M.R_V.gettrace()[-n_samples_analysis:][:,n].std(),
            ebv_values[n],
            ebv_errs[n],
            M.EBV.gettrace()[-n_samples_analysis:][:,n].mean(),
            M.EBV.gettrace()[-n_samples_analysis:][:,n].std(),
            prior_mu_info["mu"][n],
            prior_mu_info["mu_err"][n],
            M.mus.gettrace()[-n_samples_analysis:][:,n].mean(),
            M.mus.gettrace()[-n_samples_analysis:][:,n].std() ) )
    
output_stars_file.close()


post_ebv_values = array(post_ebv_values)
post_ebv_errs = array(post_ebv_errs)
# post_RV_values = array(post_RV_values)
# post_RV_errs = array(post_RV_errs)
post_mu_values = array(post_mu_values)
post_mu_errs = array(post_mu_errs)





color_vals = linspace(1, 0, n_samples)
colors_list = []
for cv in color_vals:
    colors_list.append(make_color(cv))
    
# clf()
# n = 60
# scatter(M.R_V.gettrace()[-n_samples:][:,n], M.EBV.gettrace()[-n_samples:][:,n], alpha=0.1, color=colors_list)
# 
# m = 3
# scatter(M.M_0.gettrace()[-n_samples:][:,m], M.alpha.gettrace()[-n_samples:][:,m], alpha=0.1, color=colors_list)
# 
# n=25
# clf()
# plot(M.EBV.gettrace()[-n_samples:][:,n])
# title(all_names[n] + " E(B-V) Trace  prior = %.4f +/- %.4f" % (ebv_values[n], ebv_errs[n]))
# savefig(all_names[n] + "_EBV_trace.png")
# clf()
# plot(M.R_V.gettrace()[-n_samples:][:,n])
# title(all_names[n] + " R_V Trace  prior = %.4f +/- %.4f" % (3.1, 0.75))
# savefig(all_names[n] + "_RV_trace.png")






fig = plt.figure(figsize=(3.3, 4.0))
ax1 = fig.add_subplot(1,1,1)
# plot_title = (r"$\mu_{\rm Prior}$ vs $\mu_{\rm Posterior}$ with Residual")
# fig.suptitle(plot_title, fontsize=30)
divider = make_axes_locatable(ax1)
axHistx = divider.append_axes("top", 1.25, pad=0.1, sharex=ax1)

axHistx.plot([6.8, 12.1], [0.0, 0.0], color="black")

residual_vals = []
residual_errs = []

for n in range(len(all_names)):
    blazhko_flag = False
    if rrl_info["name"][n] in blazhko_stars:
        blazhko_flag = True
    
    rrl = [rrl_info["name"][n], prior_mu_info["mu"][n], post_mu_values[n], post_mu_errs[n], prior_mu_info["mu_err"][n], rrl_info["type"][n], blazhko_flag]



    if rrl[5] == "RRab" and not rrl[6]:
        axHistx.errorbar(rrl[1], rrl[2]-rrl[1], sqrt(rrl[3]**2 + rrl[4]**2), rrl[4],
            linestyle="none",
            color="blue",
            # markeredgecolor="blue",
            marker="s",
            markersize=5, alpha=0.7)
    if rrl[5] == "RRab" and rrl[6]:
        axHistx.errorbar(rrl[1], rrl[2]-rrl[1], sqrt(rrl[3]**2 + rrl[4]**2), rrl[4],
            linestyle="none",
            color="blue",
            # markeredgecolor="blue",
            marker="D",
            markersize=5, alpha=0.7)
    if rrl[5] == "RRc" and not rrl[6]:
        axHistx.errorbar(rrl[1], rrl[2]-rrl[1], sqrt(rrl[3]**2 + rrl[4]**2), rrl[4],
            linestyle="none",
            color="red",
            # markeredgecolor="red",
            marker="s",
            markersize=5, alpha=0.7)
    if rrl[5] == "RRc" and rrl[6]:
        axHistx.errorbar(rrl[1], rrl[2]-rrl[1], sqrt(rrl[3]**2 + rrl[4]**2), rrl[4],
            linestyle="none",
            color="red",
            # markeredgecolor="red",
            marker="D",
            markersize=5, alpha=0.7)
    residual_vals.append(rrl[2]-rrl[1])
    residual_errs.append(sqrt(rrl[3]**2 + rrl[4]**2))

for tl in axHistx.get_xticklabels():
    tl.set_visible(False)
    

ax1.plot([6.8, 12.1], [6.8, 12.1], color="black")

for n in range(len(all_names)):
    blazhko_flag = False
    if rrl_info["name"][n] in blazhko_stars:
        blazhko_flag = True
    
    rrl = [rrl_info["name"][n], prior_mu_info["mu"][n], post_mu_values[n], post_mu_errs[n], prior_mu_info["mu_err"][n], rrl_info["type"][n], blazhko_flag]

    if rrl[5] == "RRab" and not rrl[6]:
        ax1.errorbar(rrl[1], rrl[2], rrl[3], rrl[4],
            linestyle="none",
            color="blue",
            # markeredgecolor="blue",
            marker="s",
            markersize=5, alpha=0.7)
    if rrl[5] == "RRab" and rrl[6]:
        ax1.errorbar(rrl[1], rrl[2], rrl[3], rrl[4],
            linestyle="none",
            color="blue",
            # markeredgecolor="blue",
            marker="D",
            markersize=5, alpha=0.7)
    if rrl[5] == "RRc" and not rrl[6]:
        ax1.errorbar(rrl[1], rrl[2], rrl[3], rrl[4],
            linestyle="none",
            color="red",
            # markeredgecolor="red",
            marker="s",
            markersize=5, alpha=0.7)
    if rrl[5] == "RRc" and rrl[6]:
        ax1.errorbar(rrl[1], rrl[2], rrl[3], rrl[4],
            linestyle="none",
            color="red",
            # markeredgecolor="red",
            marker="D",
            markersize=5, alpha=0.7)



ax1.set_xlabel(r"$\mu_{\rm Prior}$")
ax1.set_ylabel(r"$\mu_{\rm Post}$")
axHistx.set_ylabel(r"$\mu_{\rm Post} - \mu_{\rm Prior}$")

ax1.set_xlim(6.8, 12.1)
ax1.set_ylim(6.8, 12.1)
axHistx.set_ylim(-0.7, 0.7)


# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(1)
minorLocator_x = MultipleLocator(0.2)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(1)
minorLocator_y1 = MultipleLocator(0.2)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

majorLocator_y2 = MultipleLocator(0.2)
minorLocator_y2 = MultipleLocator(0.1)
axHistx.yaxis.set_major_locator(majorLocator_y2)
axHistx.yaxis.set_minor_locator(minorLocator_y2)


# pos = [left, bottom, width, height]
axHistx.set_position([0.20, 0.89, 0.76, 0.12])
ax1.set_position(    [0.20, 0.14, 0.76, 0.82])

# ax1.set_title(plot_title)
canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "mu_mu_" + plr_db_num + "_plot.pdf", dpi=300)
close("all")




plot_mag_offsets_list = [-4.5, -3.7, -3.1,   -2.6, -2.1, -1.6, -1.7, -0.6, 0.0, 0.5, 1.0,  1.5,  2.0]

metallicity_max = rrl_info["metallicity"].max()
metallicity_min = rrl_info["metallicity"].min()
metallicity_range = metallicity_max - metallicity_min

# now create the 2D colormap, multiplying by blackness to dimm vertically
color_display = zeros((256, 256, 3))
color_display[:,:,0] = (red_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.)
color_display[:,:,1] = (green_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.)
color_display[:,:,2] = (blue_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.)


# logper_grid = linspace(rrl_info["period"].min(), rrl_info["period"].max(), 1000)

logper_grid = linspace(log10(0.25/P_0), log10(0.95/P_0), 1000)

sigma_intrinsic = []
sigma_instrumental = []
sigma_prediction = []
sigma_prediction_functions = []
# in the sigma_prediction_functions, the independent variable is in logper_grid
# to transform from logper_grid into linear period, use (10**logper_grid)*P_0

plot_mag_offsets_list = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
for m in range(len(band_list)):
    fig = plt.figure(figsize = (12, 12))
    ax1 = subplot(111)

    observed_mags = []
    observed_mag_errs = []
    band_post_mus = []
    band_post_mu_errs = []
    logpers = []
    
    for n in range(len(all_names)):
        name = rrl_info["name"][n]
        if band_list[m] in mags[name].keys():
    
            period = rrl_info["period"][n]
            logpers.append(period)
            metallicity = rrl_info["metallicity"][n]
        
            observed_apparent_mag = mags[name][band_list[m]][0]
            observed_apparent_mag_err = mags[name][band_list[m]][1]
            observed_mags.append(observed_apparent_mag)
            observed_mag_errs.append(observed_apparent_mag_err)

            blazhko_flag = False
            if rrl_info["name"][n] in blazhko_stars:
                blazhko_flag = True
            
            post_mu = post_mu_values[n]
            post_mu_err = post_mu_errs[n]
            band_post_mus.append(post_mu)
            band_post_mu_errs.append(post_mu_err)
        
            post_ebv = post_ebv_values
            post_ebv_err = post_ebv_errs
            # post_RV = post_RV_values
            # post_RV_err = post_RV_errs
        
            # M_abs = m_obs + m_obs_corr - mu
        
            # apparent_mag_correction = M.R_V.gettrace()[-n_samples_analysis:][:,n]*M.EBV.gettrace()[-n_samples_analysis:][:,n]*extinction_ccm_a[m] + M.EBV.gettrace()[-n_samples_analysis:][:,n]*extinction_ccm_b[m]
            
            apparent_mag_correction = M.EBV.gettrace()[-n_samples_analysis:][:,n] * (const_R_V*extinction_ccm_a[m] + extinction_ccm_b[m])
            
            corrected_apparent_mag = observed_apparent_mag - apparent_mag_correction
                
            abs_mag_dist = corrected_apparent_mag - M.mus.gettrace()[-n_samples_analysis:][:,n]
            abs_mag = abs_mag_dist.mean()
        
            abs_mag_err = sqrt(observed_apparent_mag_err**2 + abs_mag_dist.std()**2)
        
            linear_period = (10**period) * P_0
        
#             prior_abs_mag_uncorrected = observed_apparent_mag - prior_mu_info["mu"][n]
#             ax1.plot(log10(linear_period), prior_abs_mag_uncorrected+plot_mag_offsets_list[m], linestyle="none", marker="v", color=make_color((metallicity-metallicity_min)/metallicity_range))
#             
#             prior_abs_mag_corrected = corrected_apparent_mag.mean() - prior_mu_info["mu"][n]
#             ax1.plot(log10(linear_period), prior_abs_mag_corrected+plot_mag_offsets_list[m], linestyle="none", marker="D", color=make_color((metallicity-metallicity_min)/metallicity_range))
#             
#             ax1.plot([log10(linear_period), log10(linear_period)], [prior_abs_mag_uncorrected+plot_mag_offsets_list[m], prior_abs_mag_corrected+plot_mag_offsets_list[m]], color="blue", alpha=0.7)
#             
#             ax1.plot([log10(linear_period), log10(linear_period)], [prior_abs_mag_corrected+plot_mag_offsets_list[m], abs_mag+plot_mag_offsets_list[m]], color="red", alpha=0.7)
            
            ax1.errorbar(log10(linear_period), abs_mag+plot_mag_offsets_list[m], abs_mag_err, linestyle="none", marker="s", color=make_color((metallicity-metallicity_min)/metallicity_range))
            


    observed_mags = array(observed_mags)
    observed_mag_errs = array(observed_mag_errs)
    band_post_mus = array(band_post_mus)
    band_post_mu_errs = array(band_post_mu_errs)
    logpers = array(logpers)

    sigma_by_band_mean = M.sigma_by_band.gettrace()[-n_samples_analysis:][:,m].mean()
    sigma_by_band_std = M.sigma_by_band.gettrace()[-n_samples_analysis:][:,m].std()
    sigma_from_photometry_mean = observed_mag_errs.mean()
    sigma_from_photometry_std = observed_mag_errs.std()
    
    sigma_intrinsic.append([sigma_by_band_mean, sigma_by_band_std])
    sigma_instrumental.append([sigma_from_photometry_mean, sigma_from_photometry_std])
    
    intercept = M.M_0.gettrace()[-n_samples_analysis:][:,m].mean()     
    slope = M.alpha.gettrace()[-n_samples_analysis:][:,m].mean()

    squared_m_errs = ((observed_mags - band_post_mus) - (intercept + slope*logpers))**2
    sig_m = sqrt(squared_m_errs.mean())

    ci_err_grid = []
    for logper_val in logper_grid:
        ci_err_grid.append(sqrt(sig_m**2 + std(M.M_0.gettrace()[-n_samples_analysis:][:,m] + M.alpha.gettrace()[-n_samples_analysis:][:,m] * logper_val)**2))
    ci_err_grid = array(ci_err_grid)
    
    sigma_prediction.append(ci_err_grid.min())
    
    sigma_prediction_function = interpolate.interp1d(logper_grid, ci_err_grid, kind="linear")
    sigma_prediction_functions.append(sigma_prediction_function)
    
    best_fit_line = plot_mag_offsets_list[m] + intercept + slope*logper_grid
    
    ax1.plot(log10((10**logper_grid)*P_0), best_fit_line, color="k", linestyle="solid", label=r"$M_0=$%.4f$\pm$%.4f    $\alpha=$%.4f$\pm$%.4f" %  (intercept, M.M_0.gettrace()[-n_samples_analysis:][:,m].std(), slope, M.alpha.gettrace()[-n_samples_analysis:][:,m].std()) + "\n" + r"$\sigma_{\rm intrinsic}=$%.4f$\pm$%.4f" % (sigma_by_band_mean, sigma_by_band_std) + "\n" + r"$\sigma_{\rm instrum.}=$%.4f$\pm$%.4f"%(sigma_from_photometry_mean, sigma_from_photometry_std) + "\n" + r"${\rm min}(\sigma_{\rm pred})=$%.4f" % (ci_err_grid.min()))
    ax1.plot(log10((10**logper_grid)*P_0), best_fit_line + ci_err_grid, color="k", linestyle="--")
    ax1.plot(log10((10**logper_grid)*P_0), best_fit_line - ci_err_grid, color="k", linestyle="--")


    ax1.set_xlabel("Log(Period)")
    ax1.set_ylabel("Absolute Magnitude")
    ax1.legend(loc="upper left")
    
    full_ci_range = max(best_fit_line + ci_err_grid) - min(best_fit_line - ci_err_grid)
    
    # ax1.set_ylim(max(best_fit_line + ci_err_grid)+0.025*full_ci_range, min(best_fit_line - ci_err_grid)-0.025*full_ci_range)
    ax1.set_ylim(2, -1)
    full_p_range = log10((10**rrl_info["period"].max())*P_0) - log10((10**rrl_info["period"].min())*P_0)
    
    ax1.set_xlim(log10((10**rrl_info["period"].min())*P_0)-0.025*full_p_range, log10((10**rrl_info["period"].max())*P_0)+0.025*full_p_range)
    
#     ax2 = fig.add_subplot(2,1,2)
#     ax2.imshow(color_display, origin="lower", interpolation="lanczos", 
#         extent=[metallicity_min, 0.1*(metallicity_range), 0, 0.1], alpha=1)
#     ax2.set_xlabel("Metallicity [Fe/H]")
#     ax2.xaxis.set_label_position("bottom")
#     ax2.yaxis.set_ticks([])
    # pos = [left, bottom, width, height]
    # ax1.set_position([0.075, 0.175, 0.9, 0.75])
    # ax2.set_position([0.075, 0.025, 0.9, 0.1])
    ax1.set_title(band_list[m] + " Period--Luminosity Relation")
    canvas = FigureCanvas(fig)
    canvas.print_figure(plot_dir + band_list[m] + "_" + plr_db_num + "_plr.pdf" ,transparent=True, dpi=240)
    close("all")

sigma_intrinsic = array(sigma_intrinsic)
sigma_instrumental = array(sigma_instrumental)
sigma_prediction = array(sigma_prediction)



# UBVRIzJHK from http://www.gemini.edu/?q=node/11119
#   2MASS: http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
# hipp from Equation 1.3.23 on p58 of The Hipparcos and Tycho Catalogues, Volume 1 (Section 1.3)
#       hipp catalog here: https://www.rssd.esa.int/SA/HIPPARCOS/docs/vol1_all.pdf
# W1W2W3 from http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/figures/sec4_3gt4.gif
# Bessel review article: Bessel, M. S. 2005, ARA&A, 43, 293
zero_mag_flux_dict = {  "U":1895.800,
                        "B":4266.700,
                        "hipp":3836.300, # Same as V-band, I think.
                        "V":3836.300,
                        "R":2841.000,
                        "I":2247.400,
                        "z":3564.727,
                        "J":1594.000,
                        "H":1024.000,
                        "K":666.700,
                        "W1":309.450,
                        "W2":171.787,
                        "W3":31.674
                    }



# 1 jansky = 10**(-23) erg/(s * cm**2 * Hz)
# at R=10 pc, sphere surface area is 1.19649725233e40 cm**2
jansky_lum_conversion_factor = 1.19649725233e40

L_sol = 3.839e33 # erg/s


# c = lambda * nu
# nu = c / lambda
# frequency of lambda = 1 micron is 2.99792458e14

flux_to_specral_luminosity_factor = (10**-23) * jansky_lum_conversion_factor / L_sol * 2.99792458e14


fig = plt.figure(figsize = (12, 12))
ax1 = subplot(111)
sed_period_list = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
for sed_period in sed_period_list:
    sed_mags = []
    sed_mag_errs = []
    sed_mag_errs_naive = []
    sed_wavelengths = []
    for m in range(len(plr_fitted_values)):
        band_soln = plr_fitted_values[m]
        band_name = band_soln[0]
        band_wavelength = band_soln[1]
        band_intercept = band_soln[2]
        band_intercept_err = band_soln[3]
        band_slope = band_soln[4]
        band_slope_err = band_soln[5]
        
        prediction_error_in_mags = sigma_prediction_functions[m](log10(sed_period/P_0)).mean()
        
        band_mag_at_sed_period = band_intercept + band_slope*log10(sed_period/P_0)
        band_mag_at_sed_period_err = sqrt(band_intercept_err**2 + (band_slope_err*log10(sed_period/P_0))**2)

        
        band_flux_at_sed_period = zero_mag_flux_dict[band_name] * 10**(band_mag_at_sed_period/-2.5)
    
        bright_flux_bound = zero_mag_flux_dict[band_name] * 10**((band_mag_at_sed_period-band_mag_at_sed_period_err)/-2.5)
        dim_flux_bound = zero_mag_flux_dict[band_name] * 10**((band_mag_at_sed_period+band_mag_at_sed_period_err)/-2.5)
        band_flux_at_sed_period_err = 0.5*((band_flux_at_sed_period-dim_flux_bound) + (bright_flux_bound-band_flux_at_sed_period))
        
        band_ab_mag_at_sed_period = -2.5*log10(band_flux_at_sed_period/3631.)
        
        bright_ab_mag_bound = -2.5*log10((band_flux_at_sed_period+band_flux_at_sed_period_err)/3631.)
        dim_ab_mag_bound = -2.5*log10((band_flux_at_sed_period-band_flux_at_sed_period_err)/3631.)
        band_ab_mag_at_sed_period_err = 0.5*((band_ab_mag_at_sed_period-bright_ab_mag_bound) + (dim_ab_mag_bound-band_ab_mag_at_sed_period))
        
        # sed_fluxes.append(band_flux_at_sed_period*flux_to_specral_luminosity_factor)
        # sed_flux_errs.append(band_flux_at_sed_period_err*flux_to_specral_luminosity_factor)

        sed_mags.append(band_ab_mag_at_sed_period)
        
        sed_mag_errs_naive.append(band_ab_mag_at_sed_period_err)
        sed_mag_errs.append(prediction_error_in_mags)
        
        sed_wavelengths.append(band_wavelength)

    # midline = ax1.errorbar(log10(array(sed_wavelengths)), sed_mags, sed_mag_errs, label="P=%.1f" % sed_period)
    midline = ax1.plot(log10(array(sed_wavelengths)), sed_mags, label="P=%.1f" % sed_period)
    ax1.fill_between(log10(array(sed_wavelengths)), array(sed_mags)+array(sed_mag_errs), array(sed_mags)-array(sed_mag_errs), alpha=0.25, color=getp(midline[0], "color"))

ax1.legend(loc="upper right", numpoints=1)
# ax1.set_xlim(2.8, 0)
ax1.set_xlim(log10(array(sed_wavelengths)).min()-0.1, log10(array(sed_wavelengths)).max()+0.1)
ax1.set_ylim(5.5, 0)

# ax1.set_xlabel(r"$\lambda^{-1}$ [$\mu{\rm m}^{-1}$]")
ax1.set_xlabel(r"$\log_{10}(\lambda [\mu{\rm m}])$")
# ax1.set_ylabel(r"Flux at 10 pc [Jy]")
ax1.set_ylabel(r"Absolute Magnitude (AB)")
ax1.set_title(r"RR Lyrae Spectral Energy Distribution")
canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + plr_db_num + "_SEDs.pdf", transparent=True, dpi=240)
close("all")




# from astropy.coordinates import ICRSCoordinates, GalacticCoordinates
# from astropy import units as u

fig = plt.figure(figsize = (12, 12))
ax1 = subplot(111)
# b_list = []
# my_b_list = []
for n in range(len(all_names)):
    star_ebv_color_excess = ebv_color_excess[ebv_color_excess["name"] == all_names[n]][0]
    star_ra = star_ebv_color_excess[1]
    star_dec = star_ebv_color_excess[2]
    prior_ebv = star_ebv_color_excess[3]
    prior_ebv_err = star_ebv_color_excess[4]
    post_ebv = M.EBV.gettrace()[-n_samples_analysis:][:,n].mean()
    post_ebv_err = M.EBV.gettrace()[-n_samples_analysis:][:,n].std()
    
    
    # star_coords = ICRSCoordinates(ra=star_ra, dec=star_dec, unit=(u.degree, u.degree))
    # star_galactic_latitude = star_coords.galactic.b.degrees
    
    star_galactic_latitude = 57.2957795131*arcsin( sin(star_dec*0.01745329252)*0.460199784784 - cos(star_dec*0.01745329252)*sin( (star_ra - 282.25)*0.01745329252 )*0.887815385136 )
    
    ebv_residual = (post_ebv - prior_ebv)
    ebv_residual_error = sqrt(post_ebv_err**2 + prior_ebv_err**2)
    
    ax1.errorbar(abs(star_galactic_latitude), ebv_residual, ebv_residual_error, linestyle="none", marker="s", color="blue")


# b_list = array(b_list)
# my_b_list = array(my_b_list)

ax1.set_xlabel("Abs(Galactic Latitutde) [deg]")
ax1.set_ylabel(r"$E(B-V)_{\rm post} - E(B-V)_{\rm SF}$")
ax1.set_title(r"$E(B-V)$ Residual")
canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + plr_db_num + "_EBV_residual.pdf" ,transparent=True, dpi=240)
close("all")















sys.exit()




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




# band_name = band_list[m]
# x_data = M.M_0.gettrace()[-n_samples_analysis:][:,m]
# y_data = M.alpha.gettrace()[-n_samples_analysis:][:,m]
# x_min = plr_fitted_values[m][2]-5*plr_fitted_values[m][3]
# x_max = plr_fitted_values[m][2]+5*plr_fitted_values[m][3]
# y_min = plr_fitted_values[m][4]-5*plr_fitted_values[m][5]
# y_max = plr_fitted_values[m][4]+5*plr_fitted_values[m][5]

# Add histogram panels to the top and right, maybe use Bayesian Blocks?
def contour_plotter(plr_db_num, band_name, x_data, y_data, x_min, x_max, y_min, y_max):
    x_data = array([x_data]).T
    y_data = array([y_data]).T
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

    fig = plt.figure(figsize=(8.265, 8.0))
    ax1 = fig.add_subplot(1,1,1)
    plot_title = (band_name + r" $M_0-\alpha$ Contour Plot")
    ax1.contourf(x, y, z, levels=levels_list, origin="lower", cmap=cm.gist_yarg)
    ax1.contour(x, y, z, levels=levels_list, origin="lower", colors='k')
    # ax1.scatter(x_data, y_data, alpha=0.15, color="blue", linewidths=0)
    ax1.errorbar(x_data.mean(), y_data.mean(), y_data.std(), x_data.std(), marker="o", color="r")
    ax1.set_xlabel(r"$M_0$")
    ax1.set_ylabel(r"$\alpha$")
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    ax1.set_title(plot_title)
    canvas = FigureCanvas(fig)
    canvas.print_figure(plot_dir + band_name + "_" + plr_db_num + "_contour.png", dpi=300)
    close("all")

for m in range(len(band_list)):
    contour_plotter(plr_db_num, band_list[m],
        M.M_0.gettrace()[-n_samples_analysis:][:,m], 
        M.alpha.gettrace()[-n_samples_analysis:][:,m], 
        plr_fitted_values[m][2]-5*plr_fitted_values[m][3],
        plr_fitted_values[m][2]+5*plr_fitted_values[m][3],
        plr_fitted_values[m][4]-5*plr_fitted_values[m][5],
        plr_fitted_values[m][4]+5*plr_fitted_values[m][5])





