from scipy import stats

from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

plot_dir = "plots/"


rzcep_mu = 8.0397
rzcep_mu_err = 0.0123

# From Benedict paper
hst_mu = 8.02
hst_mu_err = 0.17

# Based on HIP07: \cite{2007ASSL..350.....V}
hipp_mu = 11.34
hipp_mu_err = 2.86

gaia_mu_err = 0.021715 # 1 per cent fractional distance error

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)


plot_width_in_rzcep_sigma = 20

mu_grid = linspace(rzcep_mu-plot_width_in_rzcep_sigma*rzcep_mu_err, rzcep_mu+plot_width_in_rzcep_sigma*rzcep_mu_err, 10000)

fig = plt.figure(figsize=(3.3, 2.5))
ax1 = subplot(111)

ax1.plot(mu_grid, stats.norm.pdf(mu_grid, loc=hipp_mu, scale=hipp_mu_err), color="r", label="Hipparcos Parallax", alpha=0.6, lw=2)
ax1.plot(mu_grid, stats.norm.pdf(mu_grid, loc=hst_mu, scale=hst_mu_err), color="b", label="HST Parallax", alpha=0.6, lw=2)
ax1.plot(mu_grid, stats.norm.pdf(mu_grid, loc=rzcep_mu, scale=rzcep_mu_err), color="g", label="This Paper", alpha=0.6, lw=2)

# ax1.plot(mu_grid, stats.norm.pdf(mu_grid, loc=rzcep_mu, scale=gaia_mu_err), color="purple", label="GAIA", alpha=0.6, lw=2)

(_, caps, _) = ax1.errorbar(rzcep_mu, stats.norm.pdf(rzcep_mu, loc=rzcep_mu, scale=rzcep_mu_err)*0.60653065971263342, 0, gaia_mu_err, color="purple", label="GAIA Error", lw=2, capsize=6)
for cap in caps:
    cap.set_markeredgewidth(2)

# ax1.legend()

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.1)
minorLocator_x = MultipleLocator(0.05)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(10)
minorLocator_y1 = MultipleLocator(5)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlim(rzcep_mu-plot_width_in_rzcep_sigma*rzcep_mu_err, rzcep_mu+plot_width_in_rzcep_sigma*rzcep_mu_err)
ax1.set_ylim(-0.5, stats.norm.pdf(rzcep_mu, loc=rzcep_mu, scale=rzcep_mu_err)+1.5)

ax1.set_xlabel(r"$\mu_{\rm RZCep}$")
ax1.set_ylabel(r"$P\left(\mu_{\rm RZCep}\right)$")


ax1.annotate("This Paper", (rzcep_mu-2*rzcep_mu_err, stats.norm.pdf(8.0225, loc=rzcep_mu, scale=rzcep_mu_err)), 
                           xytext=(rzcep_mu-5*rzcep_mu_err, 15), horizontalalignment="right",
                           arrowprops=dict(linewidth=0, facecolor='green', width=2, frac=0.22, headwidth=7, shrink=0.05))

ax1.annotate("HST\nParallax", (rzcep_mu-6*rzcep_mu_err, 2.5), 
                             xytext=(rzcep_mu-7*rzcep_mu_err, 7), horizontalalignment="right",
                             arrowprops=dict(linewidth=0, facecolor='blue', width=2, frac=0.18, headwidth=7, shrink=0.05))

ax1.annotate("Hipparcos\nParallax", (8.15, 0.3), 
                                    xytext=(8.1, 6), horizontalalignment="left",
                                    arrowprops=dict(linewidth=0, facecolor='red', width=2, frac=0.15, headwidth=7, shrink=0.05))


ax1.annotate("GAIA\nError", (rzcep_mu+0.025, stats.norm.pdf(rzcep_mu, loc=rzcep_mu, scale=rzcep_mu_err)*0.60653065971263342), 
    xytext=(rzcep_mu+0.08, stats.norm.pdf(rzcep_mu, loc=rzcep_mu, scale=rzcep_mu_err)*0.60653065971263342-2.5), 
    arrowprops=dict(linewidth=0, facecolor='purple', width=2, frac=0.2, headwidth=7, shrink=0.05))


# pos =         [left, bottom, width, height]
ax1.set_position([0.16, 0.185, 0.80, 0.80])

canvas = FigureCanvas(fig)
canvas.print_figure(plot_dir + "RZCep_mu_prediction_PDF.pdf", dpi=300)
close("all")

rzcep_fwhm = 2.35482 * rzcep_mu_err

