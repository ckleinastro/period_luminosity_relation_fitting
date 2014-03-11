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

output_file = file("full_rrl_fit_table.txt", "w")
output_file.write("name\t")
output_file.write("type\t")
output_file.write("blazhko\t")
output_file.write("period\t")
output_file.write("log10(P_f/P_0)\t")
output_file.write("metallicity\t")
output_file.write("prior_mu\t")
output_file.write("prior_mu_err\t")
output_file.write("SF_ebv\t")
output_file.write("SF_ebv_err\t")
for band in band_list:
    output_file.write("m_" + band + "_obs\t")
    output_file.write("m_" + band + "_obs_err\t")
output_file.write("post_mu\t")
output_file.write("post_mu_err\t")
output_file.write("post_ebv\t")
output_file.write("post_ebv_err\t")
for band in band_list[:-1]:
    output_file.write("M_" + band + "_fit\t")
    output_file.write("M_" + band + "_fit_err\t")
output_file.write("M_" + band_list[-1] + "_fit\t")
output_file.write("M_" + band_list[-1] + "_fit_err\n")

for m in range(len(all_names)):
    full_mu_trace = []
    full_ebv_trace = []
    for n in trace_run_nums:
        mu_trace = loadtxt(data_dir + all_names[m] + "_" + plr_db_prefix + "_1" + str(n) + "_mu_trace.txt")
        ebv_trace = loadtxt(data_dir + all_names[m] + "_" + plr_db_prefix + "_1" + str(n) + "_ebv_trace.txt")
        full_mu_trace.append(mu_trace.tolist())
        full_ebv_trace.append(ebv_trace.tolist())
    full_mu_trace = array(full_mu_trace)
    full_mu_trace = full_mu_trace.flatten()
    full_ebv_trace = array(full_ebv_trace)
    full_ebv_trace = full_ebv_trace.flatten()
    
    name = all_names[m]
    type = rrl_info["type"][m]
    blazhko_affected = name in blazhko_stars
    
    fundamentalized_log_period = rrl_info["period"][m]
    if type == "RRab":
        linear_period = 10**(fundamentalized_log_period + log10(P_0))
    if type == "RRc":
        linear_period = 10**(fundamentalized_log_period + log10(P_0) - 0.127)
        
    metallicity = rrl_info["metallicity"][m]
    
    mu_prior = prior_mu_info['mu'][m]
    mu_prior_err = prior_mu_info['mu_err'][m]
    ebv_prior = ebv_values[m]
    ebv_prior_err = ebv_errs[m]
    
    ebv_post = full_ebv_trace.mean()
    ebv_post_err = full_ebv_trace.std()
    mu_post = full_mu_trace.mean()
    mu_post_err = full_mu_trace.std()
    
    
    
    observed_apparent_mags = []     # This is a string. "---" means there was no measurement
    observed_apparent_mag_errs = []     # Same as above.
    fitted_absolute_mags = []           # Same as above.
    fitted_absolute_mag_errs = []       # Same as above.
    
    for n in range(len(band_list)):
        band = band_list[n]
        a = extinction_ccm_a[n]
        b = extinction_ccm_b[n]
        if band in mags[name].keys():
            observed_apparent_mags.append("%.4f" % mags[name][band][0])
            observed_apparent_mag_errs.append("%.4f" % mags[name][band][1])
            
            # M = m_obs - E(B-V)*(3.1*a + b) - mu
            abs_mag_dist = random.normal(mags[name][band][0], mags[name][band][1], len(full_ebv_trace)) - full_ebv_trace * (const_R_V*a + b) - full_mu_trace
            abs_mag = abs_mag_dist.mean()
            abs_mag_err = abs_mag_dist.std()
            
            fitted_absolute_mags.append("%.4f" % abs_mag)
            fitted_absolute_mag_errs.append("%.4f" % abs_mag_err)
            
        else:
            observed_apparent_mags.append("---")
            observed_apparent_mag_errs.append("---")
            fitted_absolute_mags.append("---")
            fitted_absolute_mag_errs.append("---")
    
    output_file.write(name + "\t")
    output_file.write(type + "\t")
    output_file.write(str(blazhko_affected) + "\t")
    output_file.write("%.4f\t" % linear_period)
    output_file.write("%.4f\t" % fundamentalized_log_period)
    output_file.write("%.2f\t" % metallicity)
    output_file.write("%.4f\t" % mu_prior)
    output_file.write("%.4f\t" % mu_prior_err)
    output_file.write("%.4f\t" % ebv_prior)
    output_file.write("%.4f\t" % ebv_prior_err)
    for n in range(len(band_list)):
        output_file.write(observed_apparent_mags[n] + "\t")
        output_file.write(observed_apparent_mag_errs[n] + "\t")
    output_file.write("%.4f\t" % mu_post)
    output_file.write("%.4f\t" % mu_post_err)
    output_file.write("%.4f\t" % ebv_post)
    output_file.write("%.4f\t" % ebv_post_err)
    for n in range(len(band_list[:-1])):
        output_file.write(fitted_absolute_mags[n] + "\t")
        output_file.write(fitted_absolute_mag_errs[n] + "\t")
    output_file.write(fitted_absolute_mags[len(band_list)-1] + "\t")
    output_file.write(fitted_absolute_mag_errs[len(band_list)-1] + "\n")
    
output_file.close()