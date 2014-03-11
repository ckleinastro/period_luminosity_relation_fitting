from numpy import loadtxt, linspace, array, exp, log10
from scipy import interpolate
from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib
import cPickle as pickle
import sys

from scipy import random

plr_db_prefix = "M_lowexpandedEBV"
# Removed "7" and "9" because they don't converge until well into the analysis period of the trace
trace_run_nums = [1, 2, 3, 4, 5, 6, 8] 
n_trace_runs = 7

data_dir = "/Volumes/Extra_HDD/Research/Multiband_PLR_Fitting/lowexpandedEBV_trace_data/"
band_list = ["U", "B", "hipp", "V", "R", "I", "z", "J", "H", "K", "W1", "W2", "W3"]

const_R_V = 3.1

all_names = array(['AACMi', 'ABUMa', 'AEBoo', 'AFVel', 'AFVir', 'AMTuc', 'AMVir',
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


def gaussian(x, mu, sig):
    return ((sig*2.50662827463)**-1) * exp(-1*( (x - mu)**2 ) / ( 2 * (sig**2)) )





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








rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rc('axes', labelsize=10)

output_file = file("plr_fit_parameters_table.txt", "w")
output_file.write("band\t")
output_file.write("M_0\t")
output_file.write("M_0_err\t")
output_file.write("alpha\t")
output_file.write("alpha_err\t")
output_file.write("sigma_intrinsic\t")
output_file.write("sigma_intrinsic_err\t")
output_file.write("sigma_instrumental\t")
output_file.write("sigma_instrumental_err\n")

for m in range(len(band_list)):
    band = band_list[m]
    
    full_M_0_trace = []
    full_alpha_trace = []
    full_sigma_trace = []
    for n in trace_run_nums:
        M_0_trace = loadtxt(data_dir + band + "_" + plr_db_prefix + "_1" + str(n) + "_M_0_trace.txt")
        alpha_trace = loadtxt(data_dir + band + "_" + plr_db_prefix + "_1" + str(n) + "_alpha_trace.txt")
        sigma_trace = loadtxt(data_dir + band + "_" + plr_db_prefix + "_1" + str(n) + "_sigma_trace.txt")
        full_M_0_trace.append(M_0_trace.tolist())
        full_alpha_trace.append(alpha_trace.tolist())
        full_sigma_trace.append(sigma_trace.tolist())
        
    full_M_0_trace = array(full_M_0_trace)
    full_M_0_trace = full_M_0_trace.flatten()
    full_alpha_trace = array(full_alpha_trace)
    full_alpha_trace = full_alpha_trace.flatten()
    full_sigma_trace = array(full_sigma_trace)
    full_sigma_trace = full_sigma_trace.flatten()
    
    
    sigma_instrumental = []
    for name in all_names:
        if band in mags[name].keys():
            sigma_instrumental.append(mags[name][band][1])
    sigma_instrumental = array(sigma_instrumental)
    if band == "hipp":
        output_file.write(band + "\t")
    else:
        output_file.write(band + "\t\t")
    output_file.write("%.4f\t" % full_M_0_trace.mean())
    output_file.write("%.4f\t" % full_M_0_trace.std())
    output_file.write("%.4f\t" % full_alpha_trace.mean())
    output_file.write("%.4f\t" % full_alpha_trace.std())
    output_file.write("%.4f\t" % full_sigma_trace.mean())
    output_file.write("%.4f\t" % full_sigma_trace.std())
    output_file.write("%.4f\t" % sigma_instrumental.mean())
    output_file.write("%.4f\n" % sigma_instrumental.std())
    
output_file.close()