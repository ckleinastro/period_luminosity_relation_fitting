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

plot_dir = "plots/"
R_V_parameter = 3.1
P_0 = 0.52853966619770265
fundamentalized_period = 0.413484
logper = log10(fundamentalized_period/P_0)

mags = array([   10.4516, 9.7888, 9.4727, 9.2304, 9.0369, 8.8288, 9.1478, 7.8482, 7.8634, 7.7455])
mag_errs = array([0.0052, 0.0047, 0.0075, 0.0123, 0.0090, 0.0166, 0.0260, 0.0035, 0.0028, 0.0047])
band_list =      ["U",    "B",    "hipp", "V",    "R",    "I",    "z",    "W1",   "W2",   "W3"]
sigma_upper_limits = 10*array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
band_wavelengths = {    "U":0.3663,
                        "B":0.4361,
                        "hipp":0.5170,
                        "V":0.5448,
                        "R":0.6407,
                        "I":0.7980,
                        "z":0.8896,
                        "W1":3.4,
                        "W2":4.6,
                        "W3":12.0
                    }

plr_fit_parameters_dtype=dtype([('band', str, 4),
                                ('M_0', 'float'),
                                ('M_0_err', 'float'),
                                ('alpha', 'float'),
                                ('alpha_err', 'float'),
                                ('sigma_intrinsic', 'float'),
                                ('sigma_intrinsic_err', 'float'),
                                ('sigma_instrumental', 'float'),
                                ('sigma_instrumental_err', 'float')])

plr_fit_parameters = loadtxt("plr_fit_parameters_table_no_JHK.txt", skiprows=1, dtype=plr_fit_parameters_dtype)





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

extinction_ccm_a = []
extinction_ccm_b = []
for band in band_list:
    extinction_coefficients = ccm_extinction_law(band_wavelengths[band])
    extinction_ccm_a.append(extinction_coefficients[0])
    extinction_ccm_b.append(extinction_coefficients[1])

extinction_ccm_a = array(extinction_ccm_a)
extinction_ccm_b = array(extinction_ccm_b)



mu = pymc.Uniform("mu", 0.0, 14.0, value=10.0, plot=True)



# M_0 = pymc.Normal("M_0", 
#     mu=plr_fit_parameters["M_0"], 
#     tau=1.0/((plr_fit_parameters["M_0_err"])**2), 
#     value=plr_fit_parameters["M_0"],
#     plot=False )
# 
# alpha = pymc.Normal("alpha", 
#     mu=plr_fit_parameters["alpha"], 
#     tau=1.0/((plr_fit_parameters["alpha_err"])**2), 
#     value=plr_fit_parameters["alpha"],
#     plot=False )

EBV = pymc.Uniform("EBV", 0.0, 2.0, value=1.0, plot=True)

# R_V = pymc.Uniform("RV", 0.0, 6.2, value=3.1, plot=True)

# sigma = pymc.Uniform('sigma', 0.0, 200.0, value=20, plot=True)

#model
@pymc.deterministic(plot=False)
def modelled_distance_modulus(mu=mu, EBV=EBV, extinction_ccm_a=extinction_ccm_a, extinction_ccm_b=extinction_ccm_b):
    return mu + EBV*(R_V_parameter*extinction_ccm_a + extinction_ccm_b)



#likelihood
# y = pymc.Normal('y', 
#     mu=modelled_distance_modulus, 
#     tau=1.0/(sigma**2 * ( mag_errs**2 + plr_fit_parameters["M_0_err"]**2 + (logper*plr_fit_parameters["alpha_err"])**2)), 
#     value=mags - plr_fit_parameters["M_0"] - plr_fit_parameters["alpha"]*logper, 
#     observed=True)

y = pymc.Normal('y', 
    mu=modelled_distance_modulus, 
    tau=1.0/(  plr_fit_parameters["sigma_intrinsic"]**2 + mag_errs**2 + plr_fit_parameters["M_0_err"]**2 + (logper*plr_fit_parameters["alpha_err"])**2  ), 
    value=mags - plr_fit_parameters["M_0"] - plr_fit_parameters["alpha"]*logper, 
    observed=True)




# Run the MCMC sampling
needed_vars_for_mcmc = locals()

M = pymc.MCMC(needed_vars_for_mcmc, db='pickle', dbname='RZCep_traces.pickle')

# M.use_step_method(pymc.Metropolis, M.mu)
# M.use_step_method(pymc.Metropolis, M.EBV)
# M.use_step_method(pymc.Metropolis, M.M_0)
# M.use_step_method(pymc.Metropolis, M.alpha)

# full trace run (245 secs)
mcmc_sample_n = 1100000
mcmc_sample_burn = 100000
mcmc_sample_thin = 10

# faster trace run (10x faster, 10x fewer samples)
# mcmc_sample_n = 110000
# mcmc_sample_burn = 10000
# mcmc_sample_thin = 1

M.sample(mcmc_sample_n, mcmc_sample_burn, mcmc_sample_thin)
n_samples = (mcmc_sample_n - mcmc_sample_burn)/mcmc_sample_thin # use the full trace to plot convergence
n_samples_analysis = int(n_samples / 2.0)
pymc.Matplot.plot(M)
close("all")
M.db.close()


def scatter2fract_dist(scatter_mag):
    upper = ((10**((scatter_mag+5)/5))/10)-1
    lower = 1-((10**((-scatter_mag+5)/5))/10)
    return (upper + lower)/2

print "RZCep distance modulus: %.4f +/- %.4f" % (mean(M.mu.trace[:]), std(M.mu.trace[:]))

# print "RZCep prediction distance modulus: %.4f +/- %.4f" % (mean(M.mu.trace[-n_samples_analysis:]), mean(M.sigma.trace[-n_samples_analysis:])*std(M.mu.trace[-n_samples_analysis:]))

dist = 10**(mean(M.mu.trace[:])/5. + 1)
# dist_fract_err = scatter2fract_dist(mean(M.sigma.trace[-n_samples_analysis:])*std(M.mu.trace[-n_samples_analysis:]))
dist_fract_err = scatter2fract_dist(std(M.mu.trace[:]))
dist_err = dist_fract_err * dist

print "RZCep prediction distance: %.1f +/- %.1f pc" % (dist, dist_err)


print "RZCep color excess: %.4f +/- %.4f" % (mean(M.EBV.trace[:]), std(M.EBV.trace[:]))

# print "RZCep R_V parameter: %.4f +/- %.4f" % (mean(M.R_V.trace[:]), std(M.R_V.trace[:]))



sample_size = 10000
for n in range(len(band_list)):
    mu_samples = random.normal(mags[n], mag_errs[n], sample_size) - random.normal(plr_fit_parameters[n]["M_0"], plr_fit_parameters[n]["M_0_err"], sample_size) - logper * random.normal(plr_fit_parameters[n]["alpha"], plr_fit_parameters[n]["alpha_err"], sample_size)
    print "%s\tdistance modulus: %.4f +/- %.4f" % (band_list[n], mu_samples.mean(), hypot(mu_samples.std(), plr_fit_parameters[n]['sigma_intrinsic'])), 
    print "\t extinction at E(B-V)=%.3f: %.4f" % (mean(M.EBV.trace[:]), mean(M.EBV.trace[-n_samples_analysis:])*(R_V_parameter*extinction_ccm_a[n] + extinction_ccm_b[n]))

print "Recall Benedict HST parallax gives distance modulus 8.03 +/- 0.1627."

# for n in range(len(band_list)):
#     print "%s alpha prior: %.4f +/- %.4f" % (band_list[n], plr_fit_parameters[n][3], plr_fit_parameters[n][4])
#     print "%s alpha post : %.4f +/- %.4f" % (band_list[n], mean(M.alpha.trace[:,n][-n_samples_analysis:]), std(M.alpha.trace[:,n][-n_samples_analysis:]))


# plot(M.db._traces['deviance'][-n_samples:])



sys.exit()

savetxt("RZCep_EBV_trace.txt", M.EBV.trace[:])
savetxt("RZCep_mu_trace.txt", M.mu.trace[:])


sample_size = 10000
ebv_grid = linspace(0, 0.6, 100)
distance_modulus_means = []
distance_modulus_errs = []
for ebv_val in ebv_grid:
    dist_modulus_samples_list = []
    for n in range(len(band_list)):
        dist_modulus_samples = random.normal(mags[n], mag_errs[n], sample_size) - random.normal(plr_fit_parameters["M_0"][n], plr_fit_parameters["M_0_err"][n], sample_size) - logper*random.normal(plr_fit_parameters["alpha"][n], plr_fit_parameters["alpha_err"][n], sample_size) - ebv_val*(R_V_parameter*extinction_ccm_a[n] + extinction_ccm_b[n])
        dist_modulus_samples_list.append(dist_modulus_samples)
    dist_modulus_samples_array = array(dist_modulus_samples_list).flatten()
    distance_modulus_means.append(dist_modulus_samples_array.mean())
    distance_modulus_errs.append(dist_modulus_samples_array.std())
distance_modulus_means = array(distance_modulus_means)
distance_modulus_errs = array(distance_modulus_errs)

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)
    


fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)

ax1.plot(ebv_grid, distance_modulus_means, linestyle="-")
ax1.plot(ebv_grid, distance_modulus_means+distance_modulus_errs, linestyle="--")
ax1.plot(ebv_grid, distance_modulus_means-distance_modulus_errs, linestyle="--")



# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.2)
minorLocator_x = MultipleLocator(0.1)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(1.0)
minorLocator_y1 = MultipleLocator(0.5)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlabel(r"$E(B-V)$")
ax1.set_ylabel(r"$\mu$") 

ax1.legend(fontsize=4)

# pos =         [left, bottom, width, height]
ax1.set_position([0.12, 0.195, 0.83, 0.79])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "mu_vs_ebv.pdf" , dpi=300)
close("all")

print distance_modulus_means[distance_modulus_errs==distance_modulus_errs.min()][0], distance_modulus_errs.min()
