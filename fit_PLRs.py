import pymc
from scipy import loadtxt, log10, ones, sqrt, mean, polyfit, linspace, random
from scipy import zeros, append, insert, hstack, vstack, concatenate, savetxt
from scipy import identity, dot, array
from numpy import vectorize
import sys
from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import operator
import pickle
from scipy import random

const_R_V = 3.1
make_prior_plr_plots = False
if len(sys.argv) > 1:
    plr_db_num = sys.argv[1]
else:
    plr_db_num = "test"

band_list = ["U", "B", "hipp", "V", "R", "I", "z", "J", "H", "K", "W1", "W2", "W3"]
# sigma_upper_limits = array([500, 500, 400, 350, 250, 150, 100, 100, 100, 100, 100, 100, 100])
# for the sigma being added in quadrature (not multiplied) in the likelihood
sigma_upper_limits = array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

# band_list = ["W1", "W2"]
# sigma_upper_limits = array([50, 50])



# M_0_prior_vals = array([1.089, 0.847, 0.671, 0.536, 0.359, 0.224, 0.664, -0.108, -0.314, -0.317, -0.495, -0.490, -0.537])
# alpha_prior_vals = array([-1.47, -0.97, -0.73, -1.20, -1.38, -1.97, -1.91, -1.95, -2.43, -2.68, -2.38, -2.39, -2.42])

# band_list = ["W1", "W2", "W3"]
# M_0_prior_vals = array([-0.495, -0.490, -0.537])
# alpha_prior_vals = array([-2.38, -2.39, -2.42])

# Remove ARPer, RZCep, and BNVul because they are outliers, due to lying so close to the Galactic plane
# Remove ATSer b/c only have W3 and hipparcos photometry, and W3 fit is significant outlier
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
RRab_names = array(['AACMi', 'ABUMa', 'AFVel', 'AFVir', 'AMVir', 'ANSer', 'ARHer',
       'ATAnd', 'ATVir', 'AVPeg', 'AVVir', 'AXLeo', 'BBEri', 'BCDra',
       'BHPeg', 'BKDra', 'BNPav', 'BPPav', 'BRAqr', 'BTDra', 'BVAqr',
       'CGPeg', 'CIAnd', 'DDHya', 'DNAqr', 'DXDel', 'FWLup', 'HHPup',
       'HKPup', 'IKHya', 'IOLyr', 'MSAra', 'RRCet', 'RRGem', 'RRLeo',
       'RRLyr', 'RSBoo', 'RUCet', 'RUScl', 'RVCap', 'RVCet', 'RVOct',
       'RVUMa', 'RWCnc', 'RWDra', 'RXCet', 'RXCol', 'RXEri', 'RYCol',
       'RYOct', 'RZCet', 'RZCVn', 'SAra', 'SCom', 'SSCVn', 'SSFor',
       'SSLeo', 'SSOct', 'STBoo', 'STCom', 'STLeo', 'STVir', 'SUDra',
       'SVEri', 'SVHya', 'SWAnd', 'SWAqr', 'SWDra', 'SXAqr', 'SXFor',
       'SZGem', 'TTCnc', 'TTLyn', 'TUUMa', 'TVCrB', 'TWBoo', 'TWHer',
       'TWLyn', 'TYAps', 'TZAur', 'ULep', 'UPic', 'UUCet', 'UUVir',
       'UVOct', 'UYBoo', 'UYCyg', 'UZCVn', 'V341Aql', 'V413CrA', 'V440Sgr',
       'V445Oph', 'V499Cen', 'V675Sgr', 'VInd', 'VWScl', 'VXHer', 'VXScl',
       'VYLib', 'VYSer', 'VZHer', 'WCrt', 'WCVn', 'WTuc', 'WYAnt', 'WYPav',
       'WZHya', 'XAri', 'XCrt', 'XXAnd', 'XXPup', 'XZAps', 'XZCyg',
       'XZDra', 'ZMic'], 
      dtype='|S7')


all_names = test_names
# all_names = RRab_names




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

mags = pickle.load( open( "mean_flux_mags_dict.p", "rb" ) )





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

# Build up the giant, sparse observed mags matrix and the design matrix, X, 
# which will use the fundamentalized log(P/P_0) and the extinction info
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



observed_apparent_mags = []
observed_apparent_mag_errs = []
observed_apparent_mag_band_mapping = []
observed_apparent_mag_periods = []
observed_apparent_mag_ebv = []
observed_apparent_mag_ebv_err = []
observed_apparent_mag_mu_prior = []
observed_apparent_mag_mu_prior_err = []
observed_apparent_mag_names = []

for n in range(len(band_list)):
    band = band_list[n]
    for name in all_names:
        if band in mags[name].keys():
            observed_apparent_mags.append(mags[name][band][0])
            observed_apparent_mag_errs.append(mags[name][band][1])
            observed_apparent_mag_band_mapping.append(n)
            observed_apparent_mag_periods.append(rrl_info["period"][rrl_info["name"] == name][0])
            observed_apparent_mag_ebv.append(ebv_color_excess["ebv"][ebv_color_excess["name"] == name][0])
            observed_apparent_mag_ebv_err.append(ebv_color_excess["ebv_err"][ebv_color_excess["name"] == name][0])
            observed_apparent_mag_mu_prior.append(prior_mu_info['mu'][prior_mu_info['name']==name][0])
            observed_apparent_mag_mu_prior_err.append(prior_mu_info['mu_err'][prior_mu_info['name']==name][0])
            observed_apparent_mag_names.append(name)
            
observed_apparent_mags = array(observed_apparent_mags)
observed_apparent_mag_errs = array(observed_apparent_mag_errs)
observed_apparent_mag_band_mapping = array(observed_apparent_mag_band_mapping)
observed_apparent_mag_periods = array(observed_apparent_mag_periods)
observed_apparent_mag_ebv = array(observed_apparent_mag_ebv)
observed_apparent_mag_ebv_err = array(observed_apparent_mag_ebv_err)
observed_apparent_mag_mu_prior = array(observed_apparent_mag_mu_prior)
observed_apparent_mag_mu_prior_err = array(observed_apparent_mag_mu_prior_err)
observed_apparent_mag_names = array(observed_apparent_mag_names)



if make_prior_plr_plots:

    logper_grid = linspace(rrl_info["period"].min(), rrl_info["period"].max(), 1000)

    rcParams['axes.unicode_minus'] = False
    matplotlib.rc('font', family="serif")
    matplotlib.rc('font', serif="Times New Roman")
    matplotlib.rc('xtick', labelsize=12)
    matplotlib.rc('ytick', labelsize=12)
    matplotlib.rc('axes', labelsize=18)

    for m in range(len(band_list)):
        band = band_list[m]
        apparent_mags = observed_apparent_mags[observed_apparent_mag_band_mapping==m]
        apparent_mag_errs = observed_apparent_mag_errs[observed_apparent_mag_band_mapping==m]
        prior_ebvs = observed_apparent_mag_ebv[observed_apparent_mag_band_mapping==m]
        prior_ebv_errs = observed_apparent_mag_ebv_err[observed_apparent_mag_band_mapping==m]
        prior_mus = observed_apparent_mag_mu_prior[observed_apparent_mag_band_mapping==m]
        prior_mu_errs = observed_apparent_mag_mu_prior_err[observed_apparent_mag_band_mapping==m]
        period_terms = observed_apparent_mag_periods[observed_apparent_mag_band_mapping==m]
        band_specific_names = observed_apparent_mag_names[observed_apparent_mag_band_mapping==m]
    
        fig = plt.figure(figsize = (12, 12))
        ax1 = subplot(111)
    
        # with some dust, poor color excess approximation
        dust_corrected_abs_mags = (apparent_mags - 0.85*prior_ebvs*(const_R_V*extinction_ccm_a[m] + extinction_ccm_b[m]) ) - prior_mus
        dust_corrected_abs_mag_errs = sqrt(apparent_mag_errs**2 + (0.85*prior_ebv_errs)**2 + prior_mu_errs**2)
    
        linear_period = (10**period_terms) * P_0
        ax1.errorbar(log10(linear_period), dust_corrected_abs_mags, dust_corrected_abs_mag_errs, linestyle="none", marker="s", color="blue") 
    
        polyfit_slope = []
        polyfit_intercept = []
        for iter in range(1000):
            fit_params = polyfit(period_terms, random.normal(dust_corrected_abs_mags, dust_corrected_abs_mag_errs), 1)
            polyfit_slope.append(fit_params[0])
            polyfit_intercept.append(fit_params[1])
        polyfit_slope = array(polyfit_slope)
        polyfit_intercept = array(polyfit_intercept)
    
        best_fit_line = polyfit_intercept.mean() + polyfit_slope.mean()*logper_grid

        squared_m_errs = ((dust_corrected_abs_mags) - (polyfit_intercept.mean() + polyfit_slope.mean()*period_terms))**2
        sig_m = sqrt(squared_m_errs.mean())

        ci_err_grid = []
        for logper_val in logper_grid:
            ci_err_grid.append(sqrt(sig_m**2 + std(polyfit_intercept + polyfit_slope * logper_val)**2))
        ci_err_grid = array(ci_err_grid)


        full_ci_range = max(best_fit_line + ci_err_grid) - min(best_fit_line - ci_err_grid)
    
        # ax1.set_ylim(max(best_fit_line + ci_err_grid)+0.025*full_ci_range, min(best_fit_line - ci_err_grid)-0.025*full_ci_range)
        ax1.set_ylim(2, -1)
    
        full_p_range = log10((10**logper_grid)*P_0)[-1] - log10((10**logper_grid)*P_0)[0]
    
        ax1.set_xlim(log10((10**logper_grid)*P_0)[0]-0.025*full_p_range, log10((10**logper_grid)*P_0)[-1]+0.025*full_p_range)
    
    
    
        ax1.plot(log10((10**logper_grid)*P_0), best_fit_line, color="k", linestyle="solid", label=r"$M_0=$%.4f$\pm$%.4f    $\alpha=$%.4f$\pm$%.4f" % (polyfit_intercept.mean(), polyfit_intercept.std(), polyfit_slope.mean(), polyfit_slope.std()))
        ax1.plot(log10((10**logper_grid)*P_0), best_fit_line + ci_err_grid, color="k", linestyle="--")
        ax1.plot(log10((10**logper_grid)*P_0), best_fit_line - ci_err_grid, color="k", linestyle="--")

        ax1.legend(loc="upper left")
        ax1.set_xlabel("Log(Period)")
        ax1.set_ylabel("Absolute Magnitude")
        ax1.set_title(band_list[m] + " Prior Period--Luminosity Relation")
        canvas = FigureCanvas(fig)
        canvas.print_figure(band_list[m] + "_prior_plr.pdf", bbox_inches="tight",transparent=True, dpi=240)
        close("all")





ebv_values = []
ebv_errs = []
for name in all_names:
    ebv_values.append(ebv_color_excess["ebv"][ebv_color_excess["name"]==name][0])
    ebv_errs.append(ebv_color_excess["ebv_err"][ebv_color_excess["name"]==name][0])
ebv_values = array(ebv_values)
ebv_errs = array(ebv_errs)


# I think there is something wrong with this function. It does not return the
# expected array when tested on W1 or W2
# Linked to using all 138 names or only a subset (cropping out the bad ones)
def create_sparse_sources_array(band):
    query_band = band
    query_band_names = []
    for mag_name in all_names:
        bands_with_mags = mags[mag_name].keys()
        if query_band in bands_with_mags:
            query_band_names.append(mag_name)
    names_in_band = array(sorted(query_band_names, key=lambda s: s.lower()))
    
    X_d_1 = zeros( (len(names_in_band), len(all_names)) )
    col_position=0
    row_position=0
    test_index=0
    for r in range(len(all_names)):
        if test_index<len(names_in_band) and all_names[r] == names_in_band[test_index]:
            X_d_1[row_position][col_position] = 1
            row_position += 1
            test_index += 1
        col_position += 1
    return X_d_1

X_d = create_sparse_sources_array(band_list[0])
for band in band_list[1:]:
    X_d = vstack(( X_d, create_sparse_sources_array(band)))


def num_sources_for_band(band):
    query_band = band
    query_band_names = []
    for mag_name in all_names:
        bands_with_mags = mags[mag_name].keys()
        if query_band in bands_with_mags:
            if mag_name in all_names:
                query_band_names.append(mag_name)
    names_in_band = array(sorted(query_band_names, key=lambda s: s.lower()))
    return len(names_in_band)

X_M0 = zeros(( X_d.shape[0], len(band_list) ))
prev_index = 0
for n in range(len(band_list)):
    new_index = prev_index + num_sources_for_band(band_list[n])
    X_M0[:,n][prev_index:new_index] = 1
    prev_index = new_index


def list_sources_for_band(band):
    query_band = band
    query_band_names = []
    for mag_name in all_names:
        bands_with_mags = mags[mag_name].keys()
        if query_band in bands_with_mags:
            if mag_name in all_names:
                query_band_names.append(mag_name)
    names_in_band = array(sorted(query_band_names, key=lambda s: s.lower()))
    return names_in_band

X_per = zeros(( X_d.shape[0], len(band_list) ))
prev_index = 0
for n in range(len(band_list)):
    names_in_band = list_sources_for_band(band_list[n])
    new_index = prev_index + len(names_in_band)
    band_period_list = []
    for name in names_in_band:
        band_period_list.append(rrl_info["period"][rrl_info["name"]==name][0])
    band_period_array = array(band_period_list)
    X_per[:,n][prev_index:new_index] = band_period_array
    prev_index = new_index


# X_a = create_sparse_sources_array(band_list[0]) * extinction_ccm_a[0]
# for n in range(len(band_list[1:])):
#     X_a = vstack(( X_a, extinction_ccm_a[1:][n]*create_sparse_sources_array(band_list[1:][n])))
# 
# X_b = create_sparse_sources_array(band_list[0]) * extinction_ccm_b[0]
# for n in range(len(band_list[1:])):
#     X_b = vstack(( X_b, extinction_ccm_b[1:][n]*create_sparse_sources_array(band_list[1:][n])))

X_combined = create_sparse_sources_array(band_list[0]) * (const_R_V*extinction_ccm_a[0] + extinction_ccm_b[0])
for n in range(len(band_list[1:])):
    X_combined = vstack(( X_combined, (const_R_V*extinction_ccm_a[1:][n] + extinction_ccm_b[1:][n])*create_sparse_sources_array(band_list[1:][n])))

# X = hstack((X_d, X_M0, X_per, X_a, X_b))

X = hstack((X_d, X_M0, X_per, X_combined))

# X = hstack((X_d, X_M0, X_per))





sys.exit()


mus = pymc.Normal("mus", 
    mu=prior_mu_info["mu"], 
    tau=1.0 / (prior_mu_info["mu_err"]**2), 
    value=prior_mu_info["mu"])



M_0 = pymc.Normal("M_0", 
    mu=zeros(len(band_list)), 
    # mu=ab_chi_sq_fit_priors,
    tau=1.0/((2.0 * ones( len(band_list) ) )**2), 
    value=zeros(len(band_list)),
    plot=True )


alpha = pymc.Normal("alpha", 
    mu=zeros(len(band_list)), 
    # mu=ab_chi_sq_fit_priors,
    tau=1.0/((5.0 * ones( len(band_list) ))**2), 
    value=zeros(len(band_list)),
    plot=True )

# R_V = pymc.Normal("R_V",
#     mu=3.1 * ones( len(prior_mu_info["mu"])),
#     tau=1.0/((0.75 * ones( len(prior_mu_info["mu"]) ))**2),
#     value=3.1 * ones( len(prior_mu_info["mu"])) )

ebv_upper_limits = []
for ebv_value in ebv_values:
    if 2.5*ebv_value > 0.125:
        ebv_upper_limits.append(2.5*ebv_value)
    else:
        ebv_upper_limits.append(0.125)

if plr_db_num.split("_")[1] == "simpleEBV":
                                                  # make this 2.5x
    EBV = pymc.Uniform("EBV", zeros(len(ebv_values)), 2.5*ebv_values, value=ebv_values )
    
if plr_db_num.split("_")[1] == "lowexpandedEBV":
                                                  # make this 2.5x
    EBV = pymc.Uniform("EBV", zeros(len(ebv_values)), ebv_upper_limits, value=ebv_values )
    
# make 13 sigmas, one for each band
# sigma = pymc.Uniform("sigma", 0.001, 50, value=5)
# Also try with a single sigma

sigma_by_band = pymc.Uniform("sigma_by_band", 0.00001*ones(len(band_list)), sigma_upper_limits, value=0.25*sigma_upper_limits)

# b = hstack((prior_mu_info["mu"], zeros(len(band_list)), zeros(len(band_list)), 3.1 * ones( len(prior_mu_info["mu"]))*ebv_values, ebv_values))
# b = hstack((prior_mu_info["mu"], zeros(len(band_list)), zeros(len(band_list)), ebv_values))
# b = hstack(( prior_mu_info["mu"], M_0_prior_vals, alpha_prior_vals ))

#model
@pymc.deterministic(plot=False)
def modelled_apparent_mags(X=X, mus=mus, M_0=M_0, alpha=alpha, EBV=EBV):
    return dot(X, hstack((mus, M_0, alpha, EBV)))
# @pymc.deterministic(plot=False)
# def modelled_apparent_mags(X=X, mus=mus, M_0=M_0, alpha=alpha, R_V=R_V, EBV=EBV):
#     return dot(X, hstack((mus, M_0, alpha, R_V*EBV, EBV)))
# @pymc.deterministic(plot=False)
# def modelled_apparent_mags(X=X, mus=mus, M_0=M_0, alpha=alpha):
#     return dot(X, hstack((mus, M_0, alpha)))

def broadcast_sigma(sigma_by_band):
    new_sigma=[]
    for k in observed_apparent_mag_band_mapping: new_sigma.append(sigma_by_band[k])
    return array(new_sigma)


#likelihood
# y = pymc.Normal('y', 
#     mu=modelled_apparent_mags, 
#     tau=1.0/( (broadcast_sigma(sigma_by_band) * observed_apparent_mag_errs)**2 ), 
#     value=observed_apparent_mags, 
#     observed=True)
# Add the sigma_by_band in quadrature.
y = pymc.Normal('y', 
    mu=modelled_apparent_mags, 
    tau=1.0/( broadcast_sigma(sigma_by_band)**2 + observed_apparent_mag_errs**2 ), 
    value=observed_apparent_mags, 
    observed=True)
    
# y = pymc.Normal('y', 
#     mu=modelled_apparent_mags, 
#     tau=1.0/( (sigma**2) * (observed_apparent_mag_errs**2) ), 
#     value=observed_apparent_mags, 
#     observed=True)


# Run the MCMC sampling
needed_vars_for_mcmc = locals()

import socket
if socket.gethostname() == "Christopher-Kleins-MacBook-Pro.local":
    M = pymc.MCMC(needed_vars_for_mcmc, db='pickle', dbname='/Volumes/Extra_HDD/Research/Multiband_PLR_Fitting/plr_' + plr_db_num + '.pickle')
elif socket.gethostname() == "anathem":
    M = pymc.MCMC(needed_vars_for_mcmc, db='pickle', dbname='/big_data/cklein/plr_fitting_traces/plr_' + plr_db_num + '.pickle')
else:
    print "Error, not running on known host, exiting."
    sys.exit()

if plr_db_num.split("_")[0] == "AMshrink":
    M.use_step_method(pymc.AdaptiveMetropolis, M.mus, shrink_if_necessary=True)
    M.use_step_method(pymc.AdaptiveMetropolis, M.EBV, shrink_if_necessary=True)
    M.use_step_method(pymc.AdaptiveMetropolis, M.M_0, shrink_if_necessary=True)
    M.use_step_method(pymc.AdaptiveMetropolis, M.alpha, shrink_if_necessary=True)
    M.use_step_method(pymc.AdaptiveMetropolis, M.sigma_by_band, shrink_if_necessary=True)

if plr_db_num.split("_")[0] == "AMnoshrink":
    # Turn off shrink_if_necessary and see if it gets too tight in the traces
    M.use_step_method(pymc.AdaptiveMetropolis, M.mus, shrink_if_necessary=False)
    M.use_step_method(pymc.AdaptiveMetropolis, M.EBV, shrink_if_necessary=False)
    M.use_step_method(pymc.AdaptiveMetropolis, M.M_0, shrink_if_necessary=False)
    M.use_step_method(pymc.AdaptiveMetropolis, M.alpha, shrink_if_necessary=False)
    M.use_step_method(pymc.AdaptiveMetropolis, M.sigma_by_band, shrink_if_necessary=False)

# if plr_db_num.split("_")[0] == "M":
else:
    # Try explicitly using Metropolis
    M.use_step_method(pymc.Metropolis, M.mus)
    M.use_step_method(pymc.Metropolis, M.EBV)
    M.use_step_method(pymc.Metropolis, M.M_0)
    M.use_step_method(pymc.Metropolis, M.alpha)
    M.use_step_method(pymc.Metropolis, M.sigma_by_band)


# if plr_db_num == "bothshrink":
#     M.use_step_method(pymc.AdaptiveMetropolis, M.mus, shrink_if_necessary=True)
#     M.use_step_method(pymc.AdaptiveMetropolis, M.EBV, shrink_if_necessary=True)
# if plr_db_num == "mushrink":
#     M.use_step_method(pymc.AdaptiveMetropolis, M.mus, shrink_if_necessary=True)

"""
Do I need to use the shrink_if_necessary option for the EBV parameter? Previous
work shows I definitely need to use it for the mus parameter. I suspect this is
needed because both parameters have one value for each star (138), whereas the
other parameters are band-specific (13 total for each).

I have decided through some trial experiments that bothshrink is the way to go.
Nope!! Using a simple Metropolis step method is the way to go. AdaptiveMetropolis
gets stuck in precisely the wrong answer for the traces.

"""


# shrink means:
# if plr_db_num.split("_")[1] == "shrink":
#     M.use_step_method(pymc.AdaptiveMetropolis, M.EBV, shrink_if_necessary=True)
# noshrink means the line is ignored

# hopefully this takes ~30 hours
mcmc_sample_n = 350000.0*6*2
mcmc_sample_burn = 0
mcmc_sample_thin = 42

# really fast
# mcmc_sample_n = 100
# mcmc_sample_burn = 0
# mcmc_sample_thin = 1

# short run (few hours?, 2.8 hours)
# mcmc_sample_n = 350000.0
# mcmc_sample_burn = 0
# mcmc_sample_thin = 7

# Long run, should take 3 days, 100,000 samples saved
# mcmc_sample_n = 30000*60*14
# mcmc_sample_burn = 0
# mcmc_sample_thin = 3*6*14

# Long run, should take 6 days, 200,000 samples saved
# mcmc_sample_n = 30000*60*14*2
# mcmc_sample_burn = 0
# mcmc_sample_thin = 3*6*14


M.sample(mcmc_sample_n, mcmc_sample_burn, mcmc_sample_thin)
n_samples = (mcmc_sample_n - mcmc_sample_burn)/mcmc_sample_thin # use the full trace to plot convergence
n_samples_analysis = int(n_samples / 2.0)
# pymc.Matplot.plot(M)
close("all")
M.db.close()




print("#band\tM_0\tM_0_err\t\talpha\talpha_err\tsigma\tsigma_err")

for n in range(len(band_list)):
    if band_list[n] != "hipp":
        tab_str = "\t"
    else:
        tab_str = "\t"
    print("%s%s%.4f\t%.4f\t\t%.4f\t%.4f\t\t%.4f\t%.4f" % 
        (   band_list[n], tab_str, 
            M.M_0.trace[-n_samples_analysis:][:,n].mean(),
            M.M_0.trace[-n_samples_analysis:][:,n].std(),
            M.alpha.trace[-n_samples_analysis:][:,n].mean(),
            M.alpha.trace[-n_samples_analysis:][:,n].std(),
            M.sigma_by_band.trace[-n_samples_analysis:][:,n].mean(),
            M.sigma_by_band.trace[-n_samples_analysis:][:,n].std(), ))


# sigma_post = M.sigma_by_band.trace[-n_samples:][:,10]
# 
# W1_M_0 = M.M_0.trace[-n_samples:][:,10]
# W1_alpha = M.alpha.trace[-n_samples:][:,10]
# 
# R_V_0_post = M.R_V.trace[-n_samples:][:,0]
# 
# EBV_0_post = M.EBV.trace[-n_samples:][:,0]

"""
Change sigma to be quadrature-added term in the likelihood
Treating all parameters with the Metropolis step method seems to work best

Started 18 runs on Tuesday, 3/4/14 at 8:54pm
Predicted to finish on Thursday, 3/6/14 around 9am

python fit_PLRs.py M_simpleEBV_11
python fit_PLRs.py M_simpleEBV_12
python fit_PLRs.py M_simpleEBV_13
python fit_PLRs.py M_simpleEBV_14
python fit_PLRs.py M_simpleEBV_15
python fit_PLRs.py M_simpleEBV_16
python fit_PLRs.py M_simpleEBV_17
python fit_PLRs.py M_simpleEBV_18
python fit_PLRs.py M_simpleEBV_19

python fit_PLRs.py M_lowexpandedEBV_11
python fit_PLRs.py M_lowexpandedEBV_12
python fit_PLRs.py M_lowexpandedEBV_13
python fit_PLRs.py M_lowexpandedEBV_14
python fit_PLRs.py M_lowexpandedEBV_15
python fit_PLRs.py M_lowexpandedEBV_16
python fit_PLRs.py M_lowexpandedEBV_17
python fit_PLRs.py M_lowexpandedEBV_18
python fit_PLRs.py M_lowexpandedEBV_19

"""

sys.exit()

from os import system

system("python analyze_PLRs.py " + plr_db_num)


sys.exit()




names_in_band = list_sources_for_band("W1")

mus_in_band = []
mu_errs_in_band = []
for n in range(len(prior_mu_info["mu"])):
    if prior_mu_info["name"][n] in names_in_band:
        mus_in_band.append(prior_mu_info["mu"][n])
        mu_errs_in_band.append(prior_mu_info["mu_err"][n])

plot(X_per[:129][:,0], observed_apparent_mags[:129] - mus_in_band, marker="s", linestyle="none")

# plot(M.db._traces['Metropolis_mus_adaptive_scale_factor'][:])
# plot(M.db._traces['Metropolis_sigma_adaptive_scale_factor'][:])

plot(M.db._traces['Metropolis_M_0_adaptive_scale_factor'][-n_samples:])
plot(M.db._traces['Metropolis_alpha_adaptive_scale_factor'][-n_samples:])
plot(M.db._traces['Metropolis_EBV_adaptive_scale_factor'][-n_samples:])
plot(M.db._traces['deviance'][-n_samples:])
