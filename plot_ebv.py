from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib.ticker import MaxNLocator
from numpy import array, log10, loadtxt, exp, arange, hypot
from scipy import interpolate

import sys

plot_dir = "plots/"

rzcep_b = 5.501252759250636
rzcep_l = 109.46976666024584
rzcep_sf_ebv = 0.9054
rzcep_sf_ebv_err = 0.0148

rzcep_mu = 8.0397
rzcep_mu_err = 0.0123
rzcep_ebv =  0.2461
rzcep_ebv_err = 0.0089


band_list =              ["U",  "B", "hipp",  "V",  "R",  "I",  "z",  "J", "H", "K", "W1", "W2", "W3"]
plot_mag_offsets_list = [-6.1, -5.0,   -4.3, -3.6, -2.9, -2.2, -2.1, -0.8, 0.0, 0.5,  1.1,  1.6,  2.2]

P_0 = 0.52853966619770265

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)


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


ebv_color_excess = []
# name, ra, dec, E(B-V), err_E(B-V)
sfd_extinction_file = file("extinction.tbl.txt")
for line in sfd_extinction_file:
    if line[0] not in ("|", "\\"):
        if line.split()[0] in all_names:
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



sfd_ebv_color_excess = []
# name, ra, dec, E(B-V), err_E(B-V)
sfd_extinction_file = file("extinction.tbl.txt")
for line in sfd_extinction_file:
    if line[0] not in ("|", "\\"):
        if line.split()[0] in all_names:
            sfd_ebv_color_excess.append((line.split()[0], 
                                 float(line.split()[1]), float(line.split()[2]), 
                                 float(line.split()[11]), float(line.split()[12])))
sfd_extinction_file.close()
extinction_dtype = dtype([('name', str, 7),
                          ('ra', 'float'),
                          ('dec', 'float'),
                          ('ebv', 'float'),
                          ('ebv_err', 'float')])
sfd_ebv_color_excess = array(sfd_ebv_color_excess, dtype=extinction_dtype)



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


post_names = []
fitted_ebv_values = []
fitted_ebv_errs = []
post_mus = []

star_fit_data_file = file("full_rrl_fit_table.txt", "r")
first_line = star_fit_data_file.readline()
for line in star_fit_data_file:
    post_names.append(line.split()[0])
    fitted_ebv_values.append(float(line.split()[38]))
    fitted_ebv_errs.append(float(line.split()[39]))
    post_mus.append(float(line.split()[36]))
star_fit_data_file.close()      
fitted_ebv_values = array(fitted_ebv_values)
fitted_ebv_errs = array(fitted_ebv_errs)
post_mus = array(post_mus)


from astropy.coordinates import ICRSCoordinates, GalacticCoordinates
from astropy import units as u

galactic_latitudes = []
galactic_longitudes = []

for n in range(len(all_names)):
    star_ra = ebv_color_excess[n][1]
    star_dec = ebv_color_excess[n][2]
    star_coords = ICRSCoordinates(ra=star_ra, dec=star_dec, unit=(u.degree, u.degree))
    galactic_latitudes.append(star_coords.galactic.b.degrees)
    galactic_longitudes.append(star_coords.galactic.l.degrees)

galactic_latitudes = array(galactic_latitudes)
galactic_longitudes = array(galactic_longitudes)



"""
rzcep_b = 5.501252759250636
rzcep_l = 109.46976666024584
rzcep_sf_ebv = 0.9054
rzcep_sf_ebv_err = 0.0148

rzcep_mu = 8.0397
rzcep_mu_err = 0.0123
rzcep_ebv =  0.2461
rzcep_ebv_err = 0.0089
"""




fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)

ax1.hlines(0, 0, 90, linestyle="dashed", color="gray")

ebv_residual = fitted_ebv_values - ebv_color_excess["ebv"]

ebv_residual_error = sqrt(fitted_ebv_errs**2 + ebv_color_excess["ebv_err"]**2)

ax1.errorbar(abs(galactic_latitudes), ebv_residual, ebv_residual_error, linestyle="none", marker="o", ms=3, color="k")

ax1.errorbar(rzcep_b, rzcep_ebv-rzcep_sf_ebv, hypot(rzcep_ebv_err, rzcep_sf_ebv_err), linestyle="none", marker="o", ms=3, color="green")

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(10)
minorLocator_x = MultipleLocator(5)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.2)
minorLocator_y1 = MultipleLocator(0.1)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)



ax1.set_xlabel(r"$|b|$ [deg]")
ax1.set_ylabel(r"$E(B-V)_{\rm Post} - E(B-V)_{\rm SF}$") 

ax1.set_ylim(-0.7, 0.175)

mean_ebv_sigma = (1.0/sum((1.0/ebv_residual_error[abs(galactic_latitudes)>30])**2.0))**0.5


text_string = r"At $b>30$ deg, the posterior"
text_string_2 = r"$E(B-V)$ is $%.3f$ ($\pm %.3f$)" % (mean(ebv_residual[abs(galactic_latitudes)>30]), mean_ebv_sigma)
text_string_3 = r"larger than the prior SF value." 
vspacing = 0.07
ax1.text(25, -0.14, text_string, fontsize=10, ha='left', va='top')
ax1.text(25, -0.14-vspacing, text_string_2, fontsize=10, ha='left', va='top')
ax1.text(25, -0.14-2*vspacing, text_string_3, fontsize=10, ha='left', va='top')

ax1.annotate("RZCep", (rzcep_b+2.0, rzcep_ebv-rzcep_sf_ebv + 0.01), fontsize=10,
                             xytext=(rzcep_b+6, (rzcep_ebv-rzcep_sf_ebv)+0.06), horizontalalignment="left",
                             arrowprops=dict(linewidth=0, facecolor='green', width=2, frac=0.3, headwidth=7, shrink=0.05))

# pos =         [left, bottom, width, height]
ax1.set_position([0.19, 0.195, 0.77, 0.79])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "EBV_residual.pdf" , dpi=300)
close("all")


print "At b>30 deg, the fitted E(B-V) is %.4f (+/-%.4f) larger than the prior SF value, with scatter about the mean of %.4f." % (mean(ebv_residual[abs(galactic_latitudes)>30]), mean_ebv_sigma, std(ebv_residual[abs(galactic_latitudes)>30]))

ebv_residual_error[abs(galactic_latitudes)>30]

sys.exit()


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
extinction_ccm_a = array(extinction_ccm_a)
extinction_ccm_b = array(extinction_ccm_b)

median_prior_A_lambdas = []
median_fitted_A_lambdas = []
extinction_residuals = []
median_extinction_residuals = []
for m in range(len(band_list)):
    prior_A_lambdas = ebv_color_excess["ebv"][abs(galactic_latitudes)>30] * (3.1*extinction_ccm_a[m] + extinction_ccm_b[m])
    fitted_A_lambdas = fitted_ebv_values[abs(galactic_latitudes)>30] * (3.1*extinction_ccm_a[m] + extinction_ccm_b[m])
    extinction_residuals.append(fitted_A_lambdas-prior_A_lambdas)
    median_prior_A_lambdas.append(median(prior_A_lambdas))
    median_fitted_A_lambdas.append(median(fitted_A_lambdas))
    median_extinction_residuals.append(median(fitted_A_lambdas-prior_A_lambdas))









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

# now create the 2D colormap, multiplying by blackness to dimm vertically
color_display = zeros((256, 256, 3))
color_display[:,:,0] = (red_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.)
color_display[:,:,1] = (green_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.)
color_display[:,:,2] = (blue_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.)




fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(121)
ax2 = subplot(122)

for idx in range(len(post_mus)):
    # ax1.errorbar(abs(galactic_latitudes)[idx], ebv_residual[idx], ebv_residual_error[idx], linestyle="none", marker="o", ms=3, color=make_color((post_mus[idx]-post_mus.min())/(post_mus.max()-post_mus.min())) )
    ax1.errorbar(abs(galactic_latitudes)[idx], ebv_residual[idx], ebv_residual_error[idx], linestyle="none", marker="o", ms=3, color=make_color((post_mus.max() - post_mus[idx])/(post_mus.max()-post_mus.min())) )

# (left, right, bottom, top)
# extent=[post_mus.min(), 0.1*(post_mus.max()-post_mus.min()), 0, 0.2]
# extent=[0, 0.2, post_mus.min(), 0.1*(post_mus.max()-post_mus.min())]
ax2.imshow(color_display.transpose((1,0,2)), origin="upper", interpolation="lanczos", extent=[0, 0.3, post_mus.min(), post_mus.max()], alpha=1)

ax2.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(10)
minorLocator_x = MultipleLocator(5)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.1)
minorLocator_y1 = MultipleLocator(0.05)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlabel(r"$|b|$ [deg]")
ax1.set_ylabel(r"$E(B-V)_{\rm Post} - E(B-V)_{\rm SF}$") 

ax2.set_ylabel(r"Distance Modulus ($\mu_{\rm Post}$)", labelpad=5) 
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")

majorLocator_y2 = MultipleLocator(1)
minorLocator_y2 = MultipleLocator(0.5)
ax2.yaxis.set_major_locator(majorLocator_y2)
ax2.yaxis.set_minor_locator(minorLocator_y2)



text_string = r"At $b>30$ deg, the posterior"
text_string_2 = r"$E(B-V)$ is $%.3f$ ($\pm %.3f$)" % (mean(ebv_residual[abs(galactic_latitudes)>30]), std(ebv_residual[abs(galactic_latitudes)>30]))
text_string_3 = r"larger than the prior SF value." 
vspacing = 0.03
ax1.text(10, -0.14, text_string, fontsize=10, ha='left', va='top')
ax1.text(10, -0.14-vspacing, text_string_2, fontsize=10, ha='left', va='top')
ax1.text(10, -0.14-2*vspacing, text_string_3, fontsize=10, ha='left', va='top')



# pos =         [left, bottom, width, height]
ax1.set_position([0.19, 0.195, 0.57, 0.79])
ax2.set_position([0.77, 0.05, 0.13, 0.9353])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "EBV_residual_distance_colored.pdf" , dpi=300)
close("all")








abs_galactic_longitudes = abs(galactic_longitudes)




fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(121)
ax2 = subplot(122)

for idx in range(len(abs_galactic_longitudes)):
    # ax1.errorbar(abs(galactic_latitudes)[idx], ebv_residual[idx], ebv_residual_error[idx], linestyle="none", marker="o", ms=3, color=make_color((post_mus[idx]-post_mus.min())/(post_mus.max()-post_mus.min())) )
    ax1.errorbar(abs(galactic_latitudes)[idx], ebv_residual[idx], ebv_residual_error[idx], linestyle="none", marker="o", ms=3, color=make_color((abs_galactic_longitudes.max() - abs_galactic_longitudes[idx])/(abs_galactic_longitudes.max()-abs_galactic_longitudes.min())) )

# (left, right, bottom, top)
# extent=[post_mus.min(), 0.1*(post_mus.max()-post_mus.min()), 0, 0.2]
# extent=[0, 0.2, post_mus.min(), 0.1*(post_mus.max()-post_mus.min())]
ax2.imshow(color_display.transpose((1,0,2)), origin="upper", interpolation="lanczos", extent=[0, 10, abs_galactic_longitudes.min(), abs_galactic_longitudes.max()], alpha=1)

ax2.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(10)
minorLocator_x = MultipleLocator(5)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.1)
minorLocator_y1 = MultipleLocator(0.05)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlabel(r"$|b|$ [deg]")
ax1.set_ylabel(r"$E(B-V)_{\rm Post} - E(B-V)_{\rm SF}$") 

ax2.set_ylabel(r"$|l|$ [deg]", labelpad=-.5) 
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")

majorLocator_y2 = MultipleLocator(30)
minorLocator_y2 = MultipleLocator(15)
ax2.yaxis.set_major_locator(majorLocator_y2)
ax2.yaxis.set_minor_locator(minorLocator_y2)



text_string = r"At $b>30$ deg, the posterior"
text_string_2 = r"$E(B-V)$ is $%.3f$ ($\pm %.3f$)" % (mean(ebv_residual[abs(galactic_latitudes)>30]), std(ebv_residual[abs(galactic_latitudes)>30]))
text_string_3 = r"larger than the prior SF value." 
vspacing = 0.03
ax1.text(10, -0.14, text_string, fontsize=10, ha='left', va='top')
ax1.text(10, -0.14-vspacing, text_string_2, fontsize=10, ha='left', va='top')
ax1.text(10, -0.14-2*vspacing, text_string_3, fontsize=10, ha='left', va='top')



# pos =         [left, bottom, width, height]
ax1.set_position([0.19, 0.195, 0.57, 0.79])
ax2.set_position([0.77, 0.05, 0.13, 0.9353])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "EBV_residual_gallong_colored.pdf" , dpi=300)
close("all")







# functions take input from 0-255 and return output from 0-255
red_function_rb = lambda x: 255*(1-160**(-x/255.)) * (1.0-(1.0/160.0))**-1.0
green_function_rb = lambda x: array(x*0.0)
blue_function_rb = lambda x: 255.0-red_function_rb(x)


fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)

xgrid = arange(255)

ax1.plot(xgrid, red_function_rb(xgrid), color="red")
ax1.plot(xgrid, green_function_rb(xgrid), color="green")
ax1.plot(xgrid, blue_function_rb(xgrid), color="blue")


ax1.set_ylim(0, 255)
ax1.set_xlim(0, 255)
ax1.set_ylabel("Color Mapped Value")
ax1.set_xlabel("Data Value")
# pos =         [left, bottom, width, height]
ax1.set_position([0.19, 0.195, 0.77, 0.78])
canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "colromap.pdf" , dpi=300)
close("all")



def make_color_rb(val):
    val *= 255
    color = [red_function_rb(val)/255., green_function_rb(val)/255., blue_function_rb(val)/255., 1.0]
    return color

# now create the 2D colormap, multiplying by blackness to dimm vertically
color_display_rb = zeros((256, 256, 3))
color_display_rb[:,:,0] = (red_function_rb(arange(256*256).reshape(256,256)/256).transpose()  / 255.)
color_display_rb[:,:,1] = (green_function_rb(arange(256*256).reshape(256,256)/256).transpose()  / 255.)
color_display_rb[:,:,2] = (blue_function_rb(arange(256*256).reshape(256,256)/256).transpose()  / 255.)







# def make_color_rb(val):
#     val *= 255
#     color = [val/255., 0.0, (255.0-val)/255., 1.0]
#     return color
# 
# # now create the 2D colormap, multiplying by blackness to dimm vertically
# color_display_rb = zeros((256, 256, 3))
# color_display_rb[:,:,0] = ((arange(256*256).reshape(256,256)/256).transpose()  / 255.)
# color_display_rb[:,:,1] = ((zeros(256*256).reshape(256,256)/256).transpose()  / 255.)
# color_display_rb[:,:,2] = ((256-arange(256*256).reshape(256,256)/256).transpose()  / 255.)



color_list = []

for idx in range(len(fitted_ebv_values)):
    color_list.append(make_color_rb( (fitted_ebv_values[idx]-fitted_ebv_values.min())/(fitted_ebv_values.max()-fitted_ebv_values.min())    ) )
color_array = array(color_list)



sizes_array = ((1.0/post_mus) - (1.0/post_mus).min()) * (100/max((1.0/post_mus) - (1.0/post_mus).min())) + 5

axsizes_mus = array([7, 8, 9, 10, 11, 12])

axsizes_array = ((1.0/axsizes_mus) - (1.0/post_mus).min()) * (100/max((1.0/post_mus) - (1.0/post_mus).min())) + 5



"""
rzcep_b = 5.501252759250636
rzcep_l = 109.46976666024584
rzcep_sf_ebv = 0.9054
rzcep_sf_ebv_err = 0.0148

rzcep_mu = 8.0397
rzcep_mu_err = 0.0123
rzcep_ebv =  0.2461
rzcep_ebv_err = 0.0089
"""

rzcep_marker_size = ((1.0/rzcep_mu) - (1.0/post_mus).min()) * (100/max((1.0/post_mus) - (1.0/post_mus).min())) + 5
rzcep_marker_color = make_color_rb( (rzcep_ebv-fitted_ebv_values.min())/(fitted_ebv_values.max()-fitted_ebv_values.min())    ) 

fig = plt.figure(figsize = (3.3, 2.7))
ax1 = subplot(311, projection = "aitoff")

ax1.grid(True)

ax1.scatter(radians(galactic_longitudes), radians(galactic_latitudes), marker="o", color=color_array, alpha=1.0, s=sizes_array, linewidth=0.5, edgecolor='black')

ax1.scatter(radians(rzcep_l), radians(rzcep_b), marker="o", color=rzcep_marker_color, alpha=1.0, s=rzcep_marker_size, linewidth=2.5, edgecolor='green')
ax1.scatter(radians(rzcep_l), radians(rzcep_b), marker="o", color=rzcep_marker_color, alpha=1.0, s=rzcep_marker_size, linewidth=0.5, edgecolor='black')


ax1.annotate("RZCep", (radians(rzcep_l+9), radians(rzcep_b)), fontsize=10,
                             xytext=(radians(rzcep_l+22), radians(rzcep_b)), ha="left", va="center",
                             arrowprops=dict(linewidth=0, facecolor='green', width=2, frac=0.3, headwidth=7, shrink=0.05))


yticks = ax1.yaxis.get_major_ticks()
for y_tick_label in yticks:
    y_tick_label.label1.set_visible(False)

xticks = ax1.xaxis.get_major_ticks()
for x_tick_label in xticks:
    x_tick_label.label1.set_visible(False)

axsizes = fig.add_subplot(2,1,2)

axsizes.scatter(axsizes_mus, ones(len(axsizes_mus)), marker="o", color="k", alpha=1.0, s=axsizes_array, linewidth=0.5, edgecolor='black')

axsizes.set_xlim(6.5, 12.5)

axsizes_yticks = axsizes.yaxis.get_major_ticks()
for y_tick_label in axsizes_yticks:
    y_tick_label.label1.set_visible(False)


axsizes.yaxis.set_ticks([])

axsizes.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are on
    top='off',         # ticks along the top edge are off
    labelbottom='on')

axsizes.set_frame_on(False)

for tick in axsizes.get_xaxis().get_major_ticks():
    tick.set_pad(-3.)
    tick.label1 = tick._get_text1()

axsizes.set_ylabel(r"$\mu$", labelpad=-6)



ax2 = fig.add_subplot(3,1,3)
ax2.imshow(color_display_rb, origin="lower", interpolation="lanczos", 
    extent=[fitted_ebv_values.min(), fitted_ebv_values.max(), 0, 0.02], alpha=1)
    
ax2.set_xlabel(r"$E(B-V)_{\rm Post}$", labelpad=3)
ax2.xaxis.set_label_position("bottom")
ax2.yaxis.set_ticks([])

ax2.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are on
    top='off',         # ticks along the top edge are off
    labelbottom='on')
majorLocator_y2 = MultipleLocator(0.05)
minorLocator_y2 = MultipleLocator(0.025)
ax2.yaxis.set_major_locator(majorLocator_y2)
ax2.yaxis.set_minor_locator(minorLocator_y2)

ytickscolorbar = ax2.yaxis.get_major_ticks()
for y_tick_label in ytickscolorbar:
    y_tick_label.label1.set_visible(False)

# pos =         [left, bottom, width, height]
ax1.set_position([0.025, 0.225, 0.95, 0.95])
axsizes.set_position([0.1, 0.29, 0.80, 0.10])
ax2.set_position([0.025, 0.125, 0.95, 0.15])
canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "EBV_residual_aitoff.pdf" , dpi=300)
close("all")




import pickle

mags = pickle.load( open( "mean_flux_mags_dict.p", "rb" ) )

ebv_values_with_color = []
ebv_errs_with_color = []
w1_minus_w2_color = []
w1_minus_w2_err = []

for m in range(len(all_names)):
    current_name = all_names[m]
    if ("W1" in mags[current_name].keys()) and ("W2" in mags[current_name].keys()):
        ebv_values_with_color.append(fitted_ebv_values[m])
        ebv_errs_with_color.append(fitted_ebv_errs[m])
        w1_minus_w2_color.append(mags[current_name]["W1"][0] - mags[current_name]["W2"][0])
        w1_minus_w2_err.append((mags[current_name]["W1"][1]**2 + mags[current_name]["W2"][1]**2)**0.5)
ebv_values_with_color = array(ebv_values_with_color)
ebv_errs_with_color = array(ebv_errs_with_color)
w1_minus_w2_color = array(w1_minus_w2_color)
w1_minus_w2_err = array(w1_minus_w2_err)




fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)

ax1.errorbar(w1_minus_w2_color, ebv_values_with_color, ebv_errs_with_color, w1_minus_w2_err, linestyle="none", marker="o", ms=3, color="k", alpha=0.5)

# ax1.scatter(w1_minus_w2_color, ebv_values_with_color, marker="o", s=3, alpha=0.5)


# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.04)
minorLocator_x = MultipleLocator(0.02)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.1)
minorLocator_y1 = MultipleLocator(0.05)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlabel(r"$W1 - W2$")
ax1.set_ylabel(r"$E(B-V)_{\rm Post}$") 


# pos =         [left, bottom, width, height]
ax1.set_position([0.17, 0.17, 0.77, 0.79])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "W1_minus_W2_vs_ebv.pdf" , dpi=300)
close("all")







ebv_values_with_color = []
ebv_errs_with_color = []
hipp_minus_w1_color = []
hipp_minus_w1_err = []
sf_ebv_values = []
sf_ebv_errs = []

hipp_minus_w1_corrected_color = []
hipp_minus_w1_corrected_err = []

for m in range(len(all_names)):
    current_name = all_names[m]
    if ("hipp" in mags[current_name].keys()) and ("W1" in mags[current_name].keys()):
        ebv_values_with_color.append(fitted_ebv_values[m])
        ebv_errs_with_color.append(fitted_ebv_errs[m])
        sf_ebv_values.append(ebv_color_excess["ebv"][m])
        sf_ebv_errs.append(ebv_color_excess["ebv_err"][m])
        hipp_minus_w1_color.append(mags[current_name]["hipp"][0] - mags[current_name]["W1"][0])
        hipp_minus_w1_err.append((mags[current_name]["hipp"][1]**2 + mags[current_name]["W1"][1]**2)**0.5)
        
        
        
ebv_values_with_color = array(ebv_values_with_color)
ebv_errs_with_color = array(ebv_errs_with_color)
hipp_minus_w1_color = array(hipp_minus_w1_color)
hipp_minus_w1_err = array(hipp_minus_w1_err)
sf_ebv_values = array(sf_ebv_values)
sf_ebv_errs = array(sf_ebv_errs)

hipp_minus_w1_corrected_color = array(hipp_minus_w1_corrected_color)
hipp_minus_w1_corrected_err = array(hipp_minus_w1_corrected_err)

fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)

ax1.errorbar(hipp_minus_w1_color, ebv_values_with_color, ebv_errs_with_color, hipp_minus_w1_err, linestyle="none", marker="o", ms=3, color="red", alpha=0.5, label=r"${\rm vs}$ $E(B-V)_{\rm Post}$")
ax1.errorbar(hipp_minus_w1_color, sf_ebv_values, sf_ebv_errs, hipp_minus_w1_err, linestyle="none", marker="o", ms=3, color="blue", alpha=0.5, label=r"${\rm vs}$ $E(B-V)_{\rm SF}$")


# ax1.scatter(hipp_minus_w1_color, ebv_values_with_color, marker="o", s=3, color="red", alpha=0.5, label=r"${\rm vs}$ $E(B-V)_{\rm Post}$")
# ax1.scatter(hipp_minus_w1_color, sf_ebv_values, marker="o", s=3, color="blue", alpha=0.5, label=r"${\rm vs}$ $E(B-V)_{\rm SF}$")

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.2)
minorLocator_x = MultipleLocator(0.1)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.1)
minorLocator_y1 = MultipleLocator(0.05)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlabel(r"$hipp - W1$ (uncorrected)")
ax1.set_ylabel(r"$E(B-V)$") 

ax1.set_ylim(0, 0.45)
ax1.set_xlim(0.8, 2.4)

ax1.legend(loc="upper left", fontsize=10, numpoints=1, handletextpad=-0.2, labelspacing=0.2, scatterpoints=1)

# pos =         [left, bottom, width, height]
ax1.set_position([0.17, 0.190, 0.77, 0.79])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "hipp_minus_W1_uncorrected_vs_ebv.pdf" , dpi=300)
close("all")



fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)

ax1.errorbar(hipp_minus_w1_color, ebv_values_with_color-sf_ebv_values, (sf_ebv_errs**2+ebv_errs_with_color**2)**0.5, hipp_minus_w1_err, linestyle="none", marker="o", ms=3, color="k", alpha=0.5)


# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.2)
minorLocator_x = MultipleLocator(0.1)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.1)
minorLocator_y1 = MultipleLocator(0.05)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlabel(r"$hipp - W1$ (uncorrected)")
ax1.set_ylabel(r"$E(B-V)_{\rm Post} - E(B-V)_{\rm SF}$") 

# ax1.set_ylim(0, 0.45)
ax1.set_xlim(0.8, 2.4)

# pos =         [left, bottom, width, height]
ax1.set_position([0.17, 0.190, 0.77, 0.79])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "hipp_minus_W1_uncorrected_vs_ebv_residual.pdf" , dpi=300)
close("all")