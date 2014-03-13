from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from numpy import array, log10, loadtxt

import sys

plot_dir = "plots/"

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

star_fit_data_file = file("full_rrl_fit_table.txt", "r")
first_line = star_fit_data_file.readline()
for line in star_fit_data_file:
    post_names.append(line.split()[0])
    fitted_ebv_values.append(float(line.split()[38]))
    fitted_ebv_errs.append(float(line.split()[39]))
star_fit_data_file.close()      
fitted_ebv_values = array(fitted_ebv_values)
fitted_ebv_errs = array(fitted_ebv_errs)
    


from astropy.coordinates import ICRSCoordinates, GalacticCoordinates
from astropy import units as u

galactic_latitudes = []

for n in range(len(all_names)):
    star_ra = ebv_color_excess[n][1]
    star_dec = ebv_color_excess[n][2]
    star_coords = ICRSCoordinates(ra=star_ra, dec=star_dec, unit=(u.degree, u.degree))
    galactic_latitudes.append(star_coords.galactic.b.degrees)

galactic_latitudes = array(galactic_latitudes)





fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)




ebv_residual = fitted_ebv_values - ebv_color_excess["ebv"]

ebv_residual_error = sqrt(fitted_ebv_errs**2 + ebv_color_excess["ebv_err"]**2)

ax1.errorbar(abs(galactic_latitudes), ebv_residual, ebv_residual_error, linestyle="none", marker="o", ms=3, color="k")

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
ax1.set_ylabel(r"$E(B-V)_{\rm post} - E(B-V)_{\rm SF}$") 

text_string = r"At $b>30$ deg, the posterior"
text_string_2 = r"$E(B-V)$ is $%.3f$ ($\pm %.3f$)" % (mean(ebv_residual[abs(galactic_latitudes)>30]), std(ebv_residual[abs(galactic_latitudes)>30]))
text_string_3 = r"larger than the prior SF value." 
vspacing = 0.03
ax1.text(25, -0.14, text_string, fontsize=10, ha='left', va='top')
ax1.text(25, -0.14-vspacing, text_string_2, fontsize=10, ha='left', va='top')
ax1.text(25, -0.14-2*vspacing, text_string_3, fontsize=10, ha='left', va='top')



# pos =         [left, bottom, width, height]
ax1.set_position([0.19, 0.195, 0.77, 0.79])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "EBV_residual.pdf" , dpi=300)
close("all")

print "At b>30 deg, the fitted E(B-V) is %.4f (+/-%.4f) larger than the prior SF value." % (mean(ebv_residual[abs(galactic_latitudes)>30]), std(ebv_residual[abs(galactic_latitudes)>30]))