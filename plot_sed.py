from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from numpy import array, log10, loadtxt
from scipy import interpolate
import pyfits

import sys

plot_dir = "plots/"

band_list =              ["U",  "B", "hipp",  "V",  "R",  "I",  "z",  "J", "H", "K", "W1", "W2", "W3"]



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

data_dir = "/Volumes/ExtraHDD/Research/Multiband_PLR_Fitting/lowexpandedEBV_trace_data/"
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


sigma_prediction_functions = []


logper_grid_bounds = array([-0.73, 0.25])
normed_logper_grid_bounds = logper_grid_bounds - log10(P_0)
logper_grid = linspace(normed_logper_grid_bounds[0], normed_logper_grid_bounds[-1], 2000)

for band in band_list:

    M_0_fit = plr_data_dict[band]["band_fit_parameters"]["M_0"]
    M_0_err_fit = plr_data_dict[band]["band_fit_parameters"]["M_0_err"]
    alpha_fit = plr_data_dict[band]["band_fit_parameters"]["alpha"]
    alpha_err_fit = plr_data_dict[band]["band_fit_parameters"]["alpha_err"]
    


# Code to fit the prediction interval
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
    
    sigma_prediction_function = interpolate.interp1d(logper_grid, ci_err_grid, kind="linear")
    sigma_prediction_functions.append(sigma_prediction_function)






def load_kurucz_models(temp='7000', r=5, gravity="g25"):
    filename="kurucz_model_data/ckm25_" + temp + ".fits"
    hdulist = pyfits.open(filename)
    hdulist[0].header
    scidata = hdulist[1].data
    
    angstromg_grid = scidata['WAVELENGTH']
    wavelength_grid = log10(scidata['WAVELENGTH'] / 10000.0)     # converted to microns
    flam = scidata[gravity]      # flux units erg/s/cm^2/A, also called FLAM
    hdulist.close()

    fnu = (3.34e4)*(angstromg_grid**2)*flam

    # To convert to observed flux at Earth, multiply by a factor of 
    # (R/D)^2 where R is the stellar radius, and D is the distance to Earth.

    d = 10 * 3.08568025e18 # 10 pc in cm
    radius = r * 69599000000 # 5 R_sol in cm

    fnu_abs = fnu * (radius/d)**2

    ab_mag = -2.5*log10(fnu_abs/3631.)

    return wavelength_grid, ab_mag

i_wavelength = 0.764
z_wavelength = 0.906

i_z_color_data = []

spline_order = 3
spline_smoothness = 0.05

fig = plt.figure(figsize=(3.3, 4.5))
ax1 = subplot(111)

wavelength_grid, model_mag = load_kurucz_models(temp="6250", r=7.0)
ax1.plot(wavelength_grid, model_mag, color="gray", lw=0.5, alpha=1)
# ax1.plot(wavelength_grid, model_mag, color="yellow", label="Model T=6250, R=7.0", alpha=0.2)

wavelength_grid, model_mag = load_kurucz_models(temp="7000", r=4.0)
ax1.plot(wavelength_grid, model_mag, color="gray", lw=0.5, alpha=1)
# ax1.plot(wavelength_grid, model_mag, color="orange", label="Model T=7000, R=4.0", alpha=0.2)

# sed_period_list = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
sed_period_list = [0.9, 0.75, 0.6, 0.45, 0.3]
# sed_period_list = linspace(0.2, 0.9, 30)
# sed_period_list = [0.53]
rrl_spectra_flux = []
rrl_spectra_angstroms = []
for sed_period in sed_period_list:
    sed_mags = []
    sed_mag_errs = []
    sed_fluxes = []
    sed_flux_errs = []
    sed_wavelengths = []
    for m in range(len(band_list)):
        band = band_list[m]
        band_wavelength = band_wavelengths[band]
        band_intercept = plr_data_dict[band]["band_fit_parameters"]["M_0"]
        band_intercept_err = plr_data_dict[band]["band_fit_parameters"]["M_0_err"]
        band_slope = plr_data_dict[band]["band_fit_parameters"]["alpha"]
        band_slope_err = plr_data_dict[band]["band_fit_parameters"]["alpha_err"]
        
        
        prediction_error_in_mags = sigma_prediction_functions[m](log10(sed_period/P_0)).mean()
        
        band_mag_at_sed_period = band_intercept + band_slope*log10(sed_period/P_0)
        band_mag_at_sed_period_err = sqrt(band_intercept_err**2 + (band_slope_err*log10(sed_period/P_0))**2)

        
        band_flux_at_sed_period = zero_mag_flux_dict[band] * 10**(band_mag_at_sed_period/-2.5)
    
        bright_flux_bound = zero_mag_flux_dict[band] * 10**((band_mag_at_sed_period-prediction_error_in_mags)/-2.5)
        dim_flux_bound = zero_mag_flux_dict[band] * 10**((band_mag_at_sed_period+prediction_error_in_mags)/-2.5)
        
        band_flux_at_sed_period_err = 0.5*((band_flux_at_sed_period-dim_flux_bound) + (bright_flux_bound-band_flux_at_sed_period))
        
        
        
        sed_fluxes.append(band_flux_at_sed_period[0])
        sed_flux_errs.append(band_flux_at_sed_period_err[0])
        
        
        band_ab_mag_at_sed_period = -2.5*log10(band_flux_at_sed_period[0]/3631.)
        sed_mags.append(band_ab_mag_at_sed_period)
        sed_mag_errs.append(prediction_error_in_mags)
        
        sed_wavelengths.append(band_wavelength)

    
    
    
    
    log_wavelength_grid = linspace(log10(array(sed_wavelengths))[0], log10(array(sed_wavelengths))[-1], 1000)
    
#     tck = interpolate.splrep(log10(array(sed_wavelengths)), sed_mags, k=3, s=0.1)
#     ax1.plot(log_wavelength_grid, interpolate.splev(log_wavelength_grid,tck,der=0), label="3, 0.1")
    
    
    
    tck = interpolate.splrep(log10(array(sed_wavelengths)), sed_mags, k=spline_order, s=spline_smoothness)
    midline = ax1.plot(log_wavelength_grid, interpolate.splev(log_wavelength_grid,tck,der=0), label="P=%.2f d" % sed_period)
    
    spline_fit_mags = interpolate.splev(log10(array(sed_wavelengths)),tck,der=0)
    
    dense_spline_fit_mags = interpolate.splev(log10(linspace(sed_wavelengths[0], sed_wavelengths[-1], 10000)),tck,der=0)
    spline_fnu = 3631. * (10**(dense_spline_fit_mags/(-2.5)))
    sed_angstroms = linspace(sed_wavelengths[0], sed_wavelengths[-1], 10000) * 10000
    rrl_spectra_flux.append(spline_fnu)
    rrl_spectra_angstroms.append(sed_angstroms)
    
    tck_up = interpolate.splrep(log10(array(sed_wavelengths)), spline_fit_mags+array(sed_mag_errs), k=spline_order, s=0.00)
    tck_down = interpolate.splrep(log10(array(sed_wavelengths)), spline_fit_mags-array(sed_mag_errs), k=spline_order, s=0.00)
        
    ax1.fill_between(log_wavelength_grid, interpolate.splev(log_wavelength_grid,tck_up,der=0), interpolate.splev(log_wavelength_grid,tck_down,der=0), alpha=0.25, color=getp(midline[0], "color"))
    
    i_mag = interpolate.splev(log10(i_wavelength),tck,der=0)
    i_mag_ul = interpolate.splev(log10(i_wavelength),tck_up,der=0)
    i_mag_ll = interpolate.splev(log10(i_wavelength),tck_down,der=0)
    i_mag_err = 0.5*((i_mag_ul - i_mag) + (i_mag - i_mag_ll))
    
    z_mag = interpolate.splev(log10(z_wavelength),tck,der=0)
    z_mag_ul = interpolate.splev(log10(z_wavelength),tck_up,der=0)
    z_mag_ll = interpolate.splev(log10(z_wavelength),tck_down,der=0)
    z_mag_err = 0.5*((z_mag_ul - z_mag) + (z_mag - z_mag_ll))
    
    i_minus_z = i_mag - z_mag
    i_minus_z_err = (i_mag_err**2 + z_mag_err**2)**0.5
    
    print "period = %.2f\t(i-z) = %.3f +/- %.3f" % (sed_period,  i_minus_z, i_minus_z_err)
    
    i_z_color_data.append((sed_period,  i_minus_z, i_minus_z_err))
    
    # midline = ax1.plot(log10(array(sed_wavelengths)), sed_mags, label="P=%.1f" % sed_period)
    # ax1.fill_between(log10(array(sed_wavelengths)), array(sed_mags)+array(sed_mag_errs), array(sed_mags)-array(sed_mag_errs), alpha=0.25, color=getp(midline[0], "color"))
    
    # midline = ax1.plot(log10(array(sed_wavelengths)), sed_fluxes, label="P=%.1f" % sed_period)
    # ax1.fill_between(log10(array(sed_wavelengths)), array(sed_fluxes)+array(sed_flux_errs), array(sed_fluxes)-array(sed_flux_errs), alpha=0.25, color=getp(midline[0], "color"))
    
    # points = ax1.errorbar(log10(array(sed_wavelengths)), sed_mags, sed_mag_errs, marker="o", linestyle="none", markersize=3, color=getp(midline[0], "color"))







ax1.legend(loc="lower left", numpoints=1, fontsize=10, handlelength=0.7, labelspacing=0.5, handletextpad=0.4)
ax1.set_xlim(log10(array(sed_wavelengths)).min()-0.1, log10(array(sed_wavelengths)).max()+0.1)
ax1.set_ylim(5.5, 0)

# ax1.set_xlabel(r"$\lambda^{-1}$ [$\mu{\rm m}^{-1}$]")
ax1.set_xlabel(r"$\log_{10}(\lambda [\mu{\rm m}])$")
# ax1.set_ylabel(r"Flux at 10 pc [Jy]")
ax1.set_ylabel(r"Absolute Magnitude (AB)")



# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.2)
minorLocator_x = MultipleLocator(0.1)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(1)
minorLocator_y1 = MultipleLocator(0.5)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)


# pos =         [left, bottom, width, height]
ax1.set_position([0.13, 0.11, 0.85, 0.87])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "seds.pdf", dpi=300)
close("all")

rrl_radii = linspace(7, 4, len(sed_period_list))
d = 10 * 3.08568025e18 # 10 pc in cm
for n in range(len(sed_period_list)):
    period = sed_period_list[n]
    wavelengths = rrl_spectra_angstroms[n]
    radius = rrl_radii[n] * 69599000000 # R_sol in cm
    fnu = rrl_spectra_flux[n] / ((radius/d)**2)
    flam = fnu / ((3.34e4)*(wavelengths**2))  # flux units erg/s/cm^2/A, also called FLAM
    # plot(wavelengths, flam)
    # plot(wavelengths, fnu)
    output_file = file("rrl_p" + str(period).split(".")[1] + "_fnu_spec.txt", "w")
    for m in range(len(wavelengths)):
        output_file.write(str(wavelengths[m]) + "\t" + str(fnu[m]) + "\n")
    output_file.close()
    

# plot(high_p_model_grid, high_p_model_flam)
# plot(low_p_model_grid, low_p_model_flam)

sys.exit()

i_z_color_data = array(i_z_color_data)
savetxt("rrl_i_z_color_data.txt", i_z_color_data)

from scipy import stats
def fit_line(x, y):
    slope, intercept, r, prob2, see = stats.linregress(x, y)
    mx = x.mean()
    sx2 = ((x-mx)**2).sum()
    sd_intercept = see * sqrt(1./len(x) + mx*mx/sx2)
    sd_slope = see * sqrt(1./sx2)
    res = (intercept + slope*x) - y
    return intercept, sd_intercept, slope, sd_slope

errorbar(i_z_color_data[:,0][6:], i_z_color_data[:,1][6:], i_z_color_data[:,2][6:])

intercept, sd_intercept, slope, sd_slope = fit_line(i_z_color_data[:,0][6:], i_z_color_data[:,1][6:])

res = (intercept + slope*i_z_color_data[:,0][6:]) - i_z_color_data[:,1][6:]
   