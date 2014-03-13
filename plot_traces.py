import pymc
from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib.ticker import NullFormatter
from numpy import array, log10, loadtxt, arange
from matplotlib.ticker import MaxNLocator

import sys

plot_dir = "plots/"

band_list =              ["U",  "B", "hipp",  "V",  "R",  "I",  "z",  "J", "H", "K", "W1", "W2", "W3"]

P_0 = 0.52853966619770265

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)


plr_db_prefix = "M_lowexpandedEBV"
# Removed "7" and "9" because they don't converge until well into the analysis period of the trace
trace_run_nums = [1, 2, 3, 4, 5, 6, 8] 



burn_in_traces_dict = {}
middle_traces_dict = {}
converged_traces_dict = {}

variable_names = [r"$\sigma_{\rm intrinsic}$", r"$M_0$", r"$\alpha$", r"$E(B-V)$"]
plot_names = ["sigma_trace", "M_0_trace", "alpha_trace", "ebv_trace"]



nullfmt = NullFormatter()


for plot_name in plot_names:
    burn_in_traces_dict[plot_name] = []
    middle_traces_dict[plot_name] = []
    converged_traces_dict[plot_name] = []


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





for t in trace_run_nums:
    plr_db_num = plr_db_prefix + "_1" + str(t)

    M = pymc.database.pickle.load('/Volumes/Extra_HDD/Research/Multiband_PLR_Fitting/plr_' + plr_db_num + '.pickle')

    n_samples = M.M_0.gettrace().shape[0]
    n_samples_analysis = int(n_samples/2.0)

    # restrict to H band
    band_num = 8
    band = band_list[band_num]

    full_sigma_trace = M.sigma_by_band.gettrace()[:][:,band_num]
    full_M_0_trace = M.M_0.gettrace()[:][:,band_num]
    full_alpha_trace = M.alpha.gettrace()[:][:,band_num]

    star_num = where(all_names=="ABUMa")[0]
    full_ebv_trace = M.EBV.gettrace()[:][:,star_num]

    burn_in_traces_dict["sigma_trace"].append(full_sigma_trace[:10000])
    middle_traces_dict["sigma_trace"].append(full_sigma_trace[10000:50000])
    converged_traces_dict["sigma_trace"].append(full_sigma_trace[50000:])

    burn_in_traces_dict["M_0_trace"].append(full_M_0_trace[:10000])
    middle_traces_dict["M_0_trace"].append(full_M_0_trace[10000:50000])
    converged_traces_dict["M_0_trace"].append(full_M_0_trace[50000:])

    burn_in_traces_dict["alpha_trace"].append(full_alpha_trace[:10000])
    middle_traces_dict["alpha_trace"].append(full_alpha_trace[10000:50000])
    converged_traces_dict["alpha_trace"].append(full_alpha_trace[50000:])

    burn_in_traces_dict["ebv_trace"].append(full_ebv_trace[:10000])
    middle_traces_dict["ebv_trace"].append(full_ebv_trace[10000:50000])
    converged_traces_dict["ebv_trace"].append(full_ebv_trace[50000:])








y_major_steps = array([0.02, 0.03, 0.4, 0.02])
x_major_steps = array([20, 10, 1, 20])


for n in range(len(plot_names)):

    fig = plt.figure(figsize=(7, 1.75))
    ax1BurnHist = fig.add_subplot(1,3,1)
    ax1Trace = fig.add_subplot(1,3,2)
    ax1ConvHist = fig.add_subplot(1,3,3)
    
    burn_in_traces = burn_in_traces_dict[plot_names[n]]
    middle_traces = middle_traces_dict[plot_names[n]]
    converged_traces = converged_traces_dict[plot_names[n]]
    
    flattened_burn_in_traces = array(burn_in_traces).flatten()
    flattened_converged_traces = array(converged_traces).flatten()
    
    hist_range = (median(flattened_converged_traces)-flattened_converged_traces.std()*4, median(flattened_converged_traces)+flattened_converged_traces.std()*4)
    
    a = ax1BurnHist.hist(flattened_burn_in_traces, bins=100, range=hist_range, normed=True, histtype="stepfilled", orientation="horizontal", color="red", alpha=1.0)
    
    for m in range(len(burn_in_traces)):
        ax1Trace.scatter(arange(10000), burn_in_traces[m], color="red", alpha=0.03, lw=0, s=1)
        ax1Trace.scatter(arange(10000, 50000), middle_traces[m], color="green", alpha=0.03, lw=0, s=1)
        ax1Trace.scatter(arange(50000, 100000), converged_traces[m], color="blue", alpha=0.03, lw=0, s=1)

    b = ax1ConvHist.hist(flattened_converged_traces, bins=100, range=hist_range, normed=True, histtype="stepfilled", orientation="horizontal", color="blue", alpha=1.0)

    ax1BurnHist.set_ylim(hist_range[0], hist_range[1])
    ax1ConvHist.set_ylim(hist_range[0], hist_range[1])
    
    hist_max = max(a[0].max(), b[0].max())
    
    ax1BurnHist.set_xlim(0, 1.1*hist_max)
    ax1ConvHist.set_xlim(0, 1.1*hist_max)

    ax1Trace.set_ylim(hist_range[0], hist_range[1])
    ax1Trace.set_xlim(0, 99999)


    ax1BurnHist.set_ylabel(variable_names[n])
    ax1BurnHist.set_xlabel("Count")
    ax1Trace.set_xlabel("Trace Iteration")
    ax1ConvHist.set_xlabel("Count")
    
    # no labels
    ax1Trace.yaxis.set_major_formatter(nullfmt)
    ax1ConvHist.yaxis.set_major_formatter(nullfmt)
    
    y_major_step = y_major_steps[n]
    x_major_step = x_major_steps[n]

    # This code draws major and minor tick lines. Major ticks get number labels.
    hist_majorLocator_y = MultipleLocator(y_major_step)
    hist_minorLocator_y = MultipleLocator(0.5*y_major_step)
    ax1BurnHist.yaxis.set_major_locator(hist_majorLocator_y)
    ax1BurnHist.yaxis.set_minor_locator(hist_minorLocator_y)
    ax1Trace.yaxis.set_major_locator(hist_majorLocator_y)
    ax1Trace.yaxis.set_minor_locator(hist_minorLocator_y)
    ax1ConvHist.yaxis.set_major_locator(hist_majorLocator_y)
    ax1ConvHist.yaxis.set_minor_locator(hist_minorLocator_y)

    hist_majorLocator_x1 = MultipleLocator(x_major_step)
    hist_minorLocator_x1 = MultipleLocator(0.5*x_major_step)
    ax1BurnHist.xaxis.set_major_locator(hist_majorLocator_x1)
    ax1BurnHist.xaxis.set_minor_locator(hist_minorLocator_x1)
    ax1ConvHist.xaxis.set_major_locator(hist_majorLocator_x1)
    ax1ConvHist.xaxis.set_minor_locator(hist_minorLocator_x1)

    trace_majorLocator_x1 = MultipleLocator(20000)
    trace_minorLocator_x1 = MultipleLocator(10000)
    ax1Trace.xaxis.set_major_locator(trace_majorLocator_x1)
    ax1Trace.xaxis.set_minor_locator(trace_minorLocator_x1)

    ax1Trace.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    ax1Trace.xaxis.set_major_locator(MaxNLocator(prune='both'))


    # pos =         [left, bottom, width, height]
    ax1BurnHist.set_position([  0.095,                      0.25,   0.145,   0.735])
    ax1Trace.set_position([     0.095+0.145+0.02,            0.25,   0.56,    0.735])
    ax1ConvHist.set_position([  0.095+0.145+0.02+0.56+0.02,  0.25,   0.145,   0.735])


    canvas = FigureCanvas(fig)
    canvas.print_figure(plot_dir + plot_names[n] + ".png" , dpi=300)
    close("all")



import cPickle as pickle
pickle.dump( burn_in_traces_dict, open( "burn_in_traces.p", "wb" ) )
pickle.dump( middle_traces_dict, open( "middle_traces.p", "wb" ) )
pickle.dump( converged_traces_dict, open( "converged_traces.p", "wb" ) )


