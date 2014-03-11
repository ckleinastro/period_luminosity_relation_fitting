import sys
from numpy import *
from os import listdir, system
from scipy import interpolate
from pylab import plt, errorbar, imshow, savefig, close
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import MultipleLocator
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from lomb_scargle_refine import lomb



object_name = "ABUMa"


plot_dir = "plots/"

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



data_directory = object_name + "_light_curves"






def bootstrap_mag(object_name, band, full_list, previous_period, n_iter=20):
    time_array = full_list[:,0] - 2440000.0
    phase_array = ((mod(time_array, previous_period))/previous_period)
    signal_array = full_list[:,1]
    signal_err_array = full_list[:,2]
    mean_flux_mag_list = []
    mean_flux_mag_err_list = []
    harmonic_models_list = []
    modeled_data_len_list = []
    for n in range(n_iter):
        try:
            resampled_signal_array = normal(signal_array, signal_err_array)
            mean_flux_mag, mean_flux_mag_err, modl_phases, modl_mags, data_len = measure_mag(time_array, resampled_signal_array, signal_err_array, previous_period)
            mean_flux_mag_list.append(mean_flux_mag)
            mean_flux_mag_err_list.append(mean_flux_mag_err)
            harmonic_models_list.append([modl_phases, modl_mags])
            modeled_data_len_list.append(data_len)
        except:
            continue
    mean_flux_mag_array = array(mean_flux_mag_list)
    mean_flux_mag_array = ma.masked_array(mean_flux_mag_array, isnan(array(mean_flux_mag_list)))
    
#     mean_flux_mag_err_array = array(mean_flux_mag_err_list)
#     mean_flux_mag_err_array = ma.masked_array(mean_flux_mag_err_array, isnan(array(mean_flux_mag_list)))
#     
#     modeled_data_len_array = array(modeled_data_len_list)
#     modeled_data_len_array = ma.masked_array(modeled_data_len_array, isnan(array(mean_flux_mag_list)))
    
    mean_flux_mag = mean_flux_mag_array.mean()
    mean_flux_mag_err = mean_flux_mag_array.std()
    
    h = array(harmonic_models_list)
    h = ma.masked_array(h, isnan(h))
    harmonic_model_std_at_phases = h[:,1,:].std(axis=0)
    harmonic_model_mean_at_phases = h[:,1,:].mean(axis=0)
    harmonic_model_mean_amp = harmonic_model_mean_at_phases.max() - harmonic_model_mean_at_phases.min()

    return phase_array, signal_array, signal_err_array, harmonic_models_list, harmonic_model_mean_at_phases

def measure_mag(time_array, signal_array, signal_err_array, previous_period):
    mag0 = signal_array.mean()
    f = 10**(-0.4*(signal_array-mag0))
    df = 0.4*log(10)*f*signal_err_array
    max_signal_model_diff=6
    while max_signal_model_diff > 5:
        flux_psd, flux_out_dict, flux_frequency_grid = pf_lomb_code(time_array, f, df, previous_period)
        mean_flux_mag = mag0-2.5*log10(flux_out_dict['trend_coef'][0])
        mean_flux_mag_err = 2.5/log(10)*flux_out_dict['trend_coef_error'][0]    
        # Clip out the +5-sigma outliers (one at a time) and refit
        signal_model_diff = abs(f - flux_out_dict['model']) / df
        max_signal_model_diff = max(signal_model_diff)
        if max_signal_model_diff > 5:
            time_array = time_array[signal_model_diff<max(signal_model_diff)]
            f = f[signal_model_diff<max(signal_model_diff)]
            df = df[signal_model_diff<max(signal_model_diff)]
    
    time_offset = flux_out_dict["time0"]
    model_phase_coefs = flux_out_dict["rel_phase"]
    model_amp_coefs = flux_out_dict["amplitude"]

    xx = previous_period*arange(1000)/999.
    modl = model_amp_coefs[0]*sin(2*pi*(xx-time_offset)/previous_period+model_phase_coefs[0])
    for i in xrange(len(model_phase_coefs)-1):
        modl += model_amp_coefs[i+1]*sin(2*pi*(xx-time_offset)*(i+2)/previous_period+model_phase_coefs[i+1])
    # modl_mags = out_dict['trend_coef'] + 2.5*log10(modl+1)
    modl_fluxes = flux_out_dict['trend_coef'] + modl
    modl_phases = xx/previous_period
    modl_mags =  mag0-2.5*log10(modl_fluxes)
    
    return mean_flux_mag, mean_flux_mag_err, modl_phases, modl_mags, len(f)


def pf_lomb_code(x, y, dy, previous_period):
    sys_err=0.05
    dy0 = sqrt(dy**2+sys_err**2)
    Xmax = x.max()
    # f0 = 1./Xmax
    # timespan
    #df = 0.1/Xmax
    df = 1/5000.00
    #df = 0.1/200 # periodogram "resolution" in frequency-space
    f0=1/previous_period - 0*df # periodogram starting (low) frequency
    #f0=df
    #fe = 10.
    fe=1/previous_period + 2*df # periodogram ending (high) frequency
    numf = int((fe-f0)/df)
    freqin = f0 + df*arange(numf,dtype='float64')
    frequency_grid = arange(f0, f0+numf*df, df)
    #frequency_grid = array([1.0/previous_period])
    #print frequency_grid
    #print 1.0/frequency_grid
# lomb(time, signal, error, f1, df, numf, nharm=8, psdmin=6., detrend_order=0,freq_zoom=10.,tone_control=1.,return_model=True,lambda0=1.,lambda0_range=[-8,6])
    psd, out_dict = lomb(x,y,dy0,f0,df,numf, nharm=4, psdmin=2000)
    return psd, out_dict, frequency_grid



filters_list =          ["U",  "B",  "hipp", "V",  "R",  "I",  "z",  "J",  "H", "K", "W1", "W2", "W3"]
plot_mag_offsets_list = [-4.1, -3.45, -2.8,   -2.4, -1.8, -1.3, -1.5, -0.6, 0.0, 0.3, 0.7,  0.9,  1.3]

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

color_vals = linspace(1, 0, len(filters_list))
colors_list = []
for cv in color_vals:
    colors_list.append(make_color(cv))







light_curve_files = []
light_curve_filelist = listdir(data_directory)
for filename in light_curve_filelist:
    if filename[-4:] == ".txt":
        light_curve_files.append(filename)


hipparcos_phase_shift = -0.24

previous_period = rrl_info["period"][rrl_info["name"]==object_name]

max_mags = []
min_mags = []

plot_title = (object_name + " Light Curves")
fig = plt.figure(figsize=(3.3, 4.5))
ax1 = fig.add_subplot(1,1,1)
# ax1.set_title(plot_title)

model_max_data = []

for f in range(len(filters_list)):
    band = filters_list[f]
    full_list = loadtxt(data_directory + "/" + object_name + "_" + band + ".txt")
    phase_array, signal_array, signal_err_array, harmonic_models_list, harmonic_model_mean_at_phases = bootstrap_mag(object_name, band, full_list, previous_period, n_iter=500)
    model_phases = harmonic_models_list[0][0]
    
    if band=="hipp":
        phase_array = phase_array + hipparcos_phase_shift
        model_phases = model_phases + hipparcos_phase_shift

    model_phases = hstack((model_phases, model_phases+1, model_phases+2, model_phases+3))
    harmonic_model_mean_at_phases = hstack((harmonic_model_mean_at_phases, harmonic_model_mean_at_phases, harmonic_model_mean_at_phases, harmonic_model_mean_at_phases))
    
    harmonic_model_mean_at_phases = harmonic_model_mean_at_phases[model_phases<=3]
    model_phases = model_phases[model_phases<=3]
    
    
    phase_array = hstack((phase_array, phase_array+1, phase_array+2, phase_array+3))
    signal_array = hstack((signal_array, signal_array, signal_array, signal_array))
    signal_err_array = hstack((signal_err_array, signal_err_array, signal_err_array, signal_err_array))
    
    signal_err_array = signal_err_array[phase_array<=3]
    signal_array = signal_array[phase_array<=3]
    phase_array = phase_array[phase_array<=3]
    
    ax1.plot(model_phases, harmonic_model_mean_at_phases+plot_mag_offsets_list[f], color=colors_list[f], linestyle="solid", lw=2)
    ax1.plot(model_phases, harmonic_model_mean_at_phases+plot_mag_offsets_list[f], color="k", linestyle="solid", lw=1)
    ax1.errorbar(phase_array, signal_array+plot_mag_offsets_list[f], signal_err_array, color=colors_list[f], marker="o", ms=3, linestyle="none")
    
    max_mags.append(max(signal_array+plot_mag_offsets_list[f] + signal_err_array))
    min_mags.append(min(signal_array+plot_mag_offsets_list[f] - signal_err_array))
        
    max_of_model = harmonic_model_mean_at_phases[where(harmonic_model_mean_at_phases==harmonic_model_mean_at_phases.min())[0][1]]+plot_mag_offsets_list[f]
    
    phase_at_model_max = model_phases[where(harmonic_model_mean_at_phases==harmonic_model_mean_at_phases.min())[0][1]]
    
    model_max_data.append([phase_at_model_max, max_of_model])
    
    text_x_position = 2.015
    text_y_position = harmonic_model_mean_at_phases[-1] + plot_mag_offsets_list[f]
    # text_string = r"$" + band + "+" + ("$%.1f" % plot_mag_offsets_list[f])
    text_string =  r"$" + band + "$"
    ax1.text(text_x_position, text_y_position, text_string, fontsize=12, ha='left', va='center')

model_max_data = array(model_max_data)

# ax1.plot(model_max_data[:,0], model_max_data[:,1], mec="k", mew=2, ls="none", marker="s", ms=10, mfc="none")

#ax1.set_xlabel("Phase (Period=%.3f day)" % previous_period)
ax1.set_xlabel("Phase")
ax1.set_ylabel("Magnitude + Offset")
mag_range = max(max_mags) - min(min_mags)
ax1.set_ylim(max(max_mags) + 0.005*mag_range, min(min_mags) - 0.005*mag_range)
ax1.set_xlim(0, 2)

majorLocator_y1 = MultipleLocator(0.5)
minorLocator_y1 = MultipleLocator(0.25)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

majorLocator_x1 = MultipleLocator(0.5)
minorLocator_x1 = MultipleLocator(0.25)
ax1.xaxis.set_major_locator(majorLocator_x1)
ax1.xaxis.set_minor_locator(minorLocator_x1)


# pos = [left, bottom, width, height]
ax1.set_position([0.185, 0.10, 0.705, 0.89])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + object_name + ".pdf", dpi=300)
close("all")


